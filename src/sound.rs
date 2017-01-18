
use num_complex::*;
use num_traits::Float;

use std::f32::consts::PI;
use std::io::{Read};
use dft;

pub use bit_vec::BitVec;

use hound;
use memory::*;

use itertools::Itertools;

use iter::{CollectChunks, CollectChunksExt};

// TODO Move to the Detector as individual field
pub const AMPLITUDE_DEVIATION_DB: f32 = 5.;
pub const PHASE_DEVIATION_DB: f32 = PI / 4.;

pub type Cplx = Complex<f32>;

use std::cmp::Ordering;
pub type Samples  = Vec<Cplx>;
pub type Spectrum = Vec<Cplx>;
pub type SpectrumSlice = [Cplx];

pub const SAMPLE_RATE:     usize = 44100;
pub const NUM_POINTS:      usize = 8192;
pub const BASE_FREQUENCY:  f32 = (SAMPLE_RATE as f32) / (NUM_POINTS as f32);

pub fn to_decibel(input: f32) -> f32 {
    20.0 * input.log10()
}

pub fn from_decibel(input: f32) -> f32 {
    (10 as f32).powf(input / 20.0)
}

#[derive(Copy, Clone, Debug, Default)]
pub struct Detector {
    pub freq:  f32, // base detector frequency
    pub band:  f32, // frequency range
    pub amp:   f32, // amplitude in dB
    pub phase: f32, // phase
}

impl Detector {
    pub fn new(freq: f32, band: f32, amp: f32, phase: f32) -> Detector {
        Detector {
            freq:  freq,
            band:  band,
            amp:   amp,
            phase: phase,
        }
    }
}

fn freq(index: usize) -> f32 {
    (index as f32) * (SAMPLE_RATE as f32) / (NUM_POINTS as f32)
}

// TODO Rewrite using https://crates.io/crates/float-cmp
fn float_cmp<F: Float>(x: F, y: F, epsilon: F) -> Ordering {
    if (x - y).abs() < epsilon { return Ordering::Equal; }

    if x < y {
        Ordering::Less
    } else {
        Ordering::Greater
    }
}

pub fn filter_detectors_inplace(spectrum: &SpectrumSlice, detectors: &[Detector], result: &mut BitVec) {
    // Iterating through all detectors filtering out activity
    for detector in detectors {
        let index = (detector.freq / BASE_FREQUENCY).round() as usize;

        println!("detector freq {} ± {}, amp {} ± {} dB, phase {}",
            detector.freq,
            detector.band,
            detector.amp,
            AMPLITUDE_DEVIATION_DB,
            detector.phase
        );

        // Selecting the entry with the largest amplitude
        let amplitude = (spectrum[index].norm() * 2.0) / NUM_POINTS as f32;
        let phase = spectrum[index].arg();

        // Treating detector as active if max amplitude lays within detector's selectivity range
        let amp_match   = (to_decibel(amplitude).abs() - detector.amp.abs()).abs() < AMPLITUDE_DEVIATION_DB;

        // + PI is required to compare positive and negative values
        let phase_match = ((phase + PI) - (detector.phase + PI)).abs() < PHASE_DEVIATION_DB;

        let is_active   = amp_match && phase_match;
        result.push(is_active);

        println!("signal amplitude {} → {} dB, phase {}{}\n",
            amplitude,
            to_decibel(amplitude),
            phase,
            if is_active { ", **MATCH**" } else { "" }
        );
    }
}

pub fn filter_detectors_inplace2(spectrum: &SpectrumSlice, detectors: &[Detector], result: &mut BitVec) {
    // Iterating through all detectors filtering out activity
    for detector in detectors {
        // Each detector operates only in the fixed part of the spectrum
        // Selecting potentially interesting spectrum slice to check
        let lo = ((detector.freq - detector.band).abs() / BASE_FREQUENCY).round() as usize;
        let hi = ((detector.freq + detector.band).abs() / BASE_FREQUENCY).round() as usize;

        if lo > spectrum.len() - 1 || hi > spectrum.len() - 1 {
            println!("invalid detector freq {}, band {}", detector.freq, detector.band);
            break;
        }

        let range = &spectrum[lo .. hi + 1];

        println!("detector freq {} ± {}, amp {} ± {} dB, phase {}, selected freq range is {} .. {}",
            detector.freq,
            detector.band,
            detector.amp,
            AMPLITUDE_DEVIATION_DB,
            detector.phase,
            freq(lo),
            freq(hi)
        );

        // Selecting the entry with the largest amplitude
        let (index, amplitude, phase) = range
            .iter()
            .enumerate()
            .map(|(i, c)| (i, (c.norm() * 2.0) / NUM_POINTS as f32, c.arg()))
            .max_by(|&(_, x, _), &(_, y, _)| float_cmp(x, y, 0.00001))
            .unwrap();

        // Treating detector as active if max amplitude lays within detector's selectivity range
        let amp_match   = (to_decibel(amplitude).abs() - detector.amp.abs()).abs() < AMPLITUDE_DEVIATION_DB;
        let phase_match = ((phase + PI) - (detector.phase + PI)).abs() < PHASE_DEVIATION_DB;
        let is_active   = amp_match && phase_match;

        println!("signal frequency {}, amplitude {} → {} dB, phase {}{}\n",
            freq(lo+index),
            amplitude,
            to_decibel(amplitude),
            phase,
            if is_active { ", **MATCH**" } else { "" }
        );

        result.push(is_active);
    }
}

pub fn analyze_file(filename: &str, detectors: &[Detector]) -> BitVec {
    // This will hold resulting bit vector of the detectors activity mask
    let mut result = BitVec::new();

    // Opening wave file for reading
    let mut reader = hound::WavReader::open(filename).unwrap();
    let mut total_samples = 0;

    // Reading file by chunks of NUM_POINTS, normalizing value to [0 .. 1], then processing
    for chunk in 1 .. {
        let mut samples: Samples = reader
            .samples::<i16>()
            .take(NUM_POINTS / 2)
            .map(|s| Cplx::new(s.unwrap() as f32 / i16::max_value() as f32, 0.0))
            .collect();

        total_samples += samples.len();

        if samples.len() < NUM_POINTS / 2 {
            println!("*** End of input, samples processed: {}, chunks {}", total_samples, chunk - 1);
            break;
        }

        println!("*** Chunk {}", chunk);
        samples.resize(NUM_POINTS, Cplx::default());

        let plan = dft::Plan::new(dft::Operation::Forward, NUM_POINTS);
        dft::transform(&mut samples, &plan);

        filter_detectors_inplace(&samples[..NUM_POINTS/2 - 1], detectors, &mut result);
        println!("total samples so far {}", total_samples);
    }

    result
}

pub fn build_dictionary<'d>(filename: &str, detectors: &'d [Detector]) -> Dictionary<'d> {
    // 1. Read 2*n samples from the stream
    // 2. Collect them into vector
    // 3. Perform series of DST on it [0 .. n] .. [n .. 2*n]

    let mut reader = hound::WavReader::open(filename).unwrap();

    let freqs = (102.28, 156.12);
    //let freqs = (900., 1100.);
    let mut dictionary = Dictionary::new(detectors, freqs.0, freqs.1);

    // Each detector operates only in the fixed part of the spectrum
    // Selecting potentially interesting spectrum slice to check
    let low = (freqs.0 / BASE_FREQUENCY).round() as usize;
    let high = (freqs.1 / BASE_FREQUENCY).round() as usize;

    let plan = dft::Plan::new(dft::Operation::Forward, NUM_POINTS);

    const SLICES_PER_FRAME:    usize = 32;
    const SLICE_OFFSET:        usize = NUM_POINTS / SLICES_PER_FRAME;
    const FRAGMENTS_PER_FRAME: usize = 4;
    const SLICES_PER_FRAGMENT: usize = (SLICES_PER_FRAME / FRAGMENTS_PER_FRAME) / 2;

    println!("Building dictionary... ");
    let mut frame_count = 0;

    let mut spectra = Vec::new();

    for frame in reader
        .samples::<i16>()
        .map(|s| Cplx::new(s.unwrap() as f32 / i16::max_value() as f32, 0.0))
        .collect_chunks::<Samples>(NUM_POINTS)
    {
        if frame.len() < NUM_POINTS {
            break;
        }

        spectra.clear();

        // N spectra spanning the whole frame shifted in time
        for slice in 0 .. SLICES_PER_FRAME / 2 {
            let range = slice * SLICE_OFFSET .. slice * SLICE_OFFSET + NUM_POINTS/2;

            let mut samples: Samples = frame[range].iter().cloned().collect();
            samples.resize(NUM_POINTS, Cplx::default());

            dft::transform(&mut samples, &plan);
            spectra.push(samples as Spectrum);
        }

        // For each registered fragment (currently the only one)
        for fragment_index in 0 .. FRAGMENTS_PER_FRAME {
            let mut fragment_spectra = Vec::new();

            // Spectrum slices for the particular fragment's frequency region
            let fragment_region = fragment_index*SLICES_PER_FRAGMENT .. (fragment_index + 1)*SLICES_PER_FRAGMENT;
            for spectrum in &spectra[fragment_region] {
                let slice: Spectrum = spectrum[low .. high + 1].iter().cloned().collect();
                fragment_spectra.push(slice);
            }

            let fragment = Fragment::from_spectra(fragment_spectra);
            dictionary.insert_fragment(fragment, 30);
        }

        frame_count += 1;
        let dictionary_size = dictionary.len();

        println!("\r{} frames processed, {} fragments classified", frame_count, dictionary_size);
    }

    println!("\nCompleted.");
    dictionary
}

pub fn generate_file(filename: &str, detectors: &[Detector], mask: &BitVec) {
    let wav_header = hound::WavSpec {
        channels: 1,
        sample_rate: SAMPLE_RATE as u32, //44100,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };

    let mut writer = hound::WavWriter::create(filename, wav_header).unwrap();
    let mut mask_iter = mask.iter();

    let mut total_samples = 0;

    for chunk in 1 .. {
        println!("!!! Chunk {}", chunk);

        let mut spectrum: Spectrum = Vec::new();
        spectrum.resize(NUM_POINTS, Cplx::default());

        let mut num_taken = 0;

        for (frequency, amplitude, phase) in mask_iter
            .clone()
            .take(detectors.len())
            .filter(|_| { num_taken += 1; true })
            .enumerate()
            .filter(|&(_, bit)| bit)
            .map(|(index, _)| {
                let frequency = detectors[index].freq;
                let amplitude = from_decibel(detectors[index].amp) / 2.0;// * NUM_POINTS as f32;
                let phase = detectors[index].phase;

                println!("chunk {}, active detector[{}]: freq {}, amplitude {} dB = {}",
                    chunk,
                    index,
                    frequency,
                    detectors[index].amp,
                    amplitude);

                (frequency, amplitude, phase)
            })
        {
            let bin = (frequency / BASE_FREQUENCY).round() as usize;

            println!("writing freq {} to bin[{}] = bin[{}] = {}\n", frequency, bin, NUM_POINTS-bin, amplitude);

            let value = Cplx::from_polar(&amplitude, &phase);

            if spectrum[bin].norm() < value.norm() {
                spectrum[bin] = value;
                // spectrum[NUM_POINTS - bin].re = amplitude;
            }
        }

        if num_taken == 0 {
            break;
        }

        // Advance iterator to match cloned position
        mask_iter.nth(num_taken - 1);

        let plan = dft::Plan::new(dft::Operation::Backward, NUM_POINTS);
        dft::transform(&mut spectrum, &plan);

        let max_sample = spectrum.iter().max_by(|&x, &y| float_cmp(x.re, y.re, 0.00001)).unwrap().re;
        //let max_sample = spectrum.iter().max_by(|x, y| x.norm().partial_cmp(&y.norm()).unwrap_or(Ordering::Less)).unwrap().re;

        for sample in spectrum {
            total_samples += 1;

            let amplitude = (i16::max_value() - 1000) as f32;
            writer.write_sample(((sample.re / max_sample) * (amplitude / 2.)) as i16).unwrap();
        }

        println!("total samples so far {}", total_samples);
    }

    println!("total samples written {}", total_samples);
}
