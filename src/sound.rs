
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

use sample::*;
use sample::window::Type;

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

const SLICES_PER_FRAME:    usize = 16;
const FRAGMENTS_PER_FRAME: usize = 4;

const SLICE_OFFSET:        usize = (NUM_POINTS / 2) / SLICES_PER_FRAME;
const SLICES_PER_FRAGMENT: usize = SLICES_PER_FRAME / FRAGMENTS_PER_FRAME;
const FRAGMENT_WINDOW: (f32, f32) = (350., 500.);

type KeyVec = Vec<Option<FragmentKey>>;

pub fn build_glossary<'d>(filename: &str, detectors: &'d [Detector]) -> (Glossary<'d>, KeyVec) {
    // (100, 199), (200, 299), ... (1900, 1999)
    let regions: Vec<_> = (1 .. 13).into_iter().map(|i: u32| (i as f32 * 100., i as f32 * 100. + 99.)).collect();
    let mut dictionaries: Vec<_> = regions.iter().map(|r| Dictionary::new(detectors, r.0, r.1)).collect();

    let plan = dft::Plan::new(dft::Operation::Forward, NUM_POINTS);
    let mut reader = hound::WavReader::open(filename).unwrap();
    let mut spectra = Vec::with_capacity(SLICES_PER_FRAME);

    println!("Building glossary... ");
    let mut frame_count = 0;

    let mut keys = KeyVec::new();

    for frame in reader
        .samples::<i16>()
        .map(|s| Cplx::new(s.unwrap() as f32 / i16::max_value() as f32, 0.0))
        .collect_chunks::<Samples>(NUM_POINTS)
    {
        if frame.len() < NUM_POINTS {
            break;
        }

        spectra.clear();

        // N spectrums spanning the whole frame shifted in time
        for slice in 0 .. SLICES_PER_FRAME {
            // Collecting samples
            let range = slice * SLICE_OFFSET .. slice * SLICE_OFFSET + NUM_POINTS / 2;
            let mut samples: Samples = frame[range].iter().cloned().collect();

            // Zero padding to the frame length
            samples.resize(NUM_POINTS, Cplx::default());

            // Performing FFT
            dft::transform(&mut samples, &plan);
            spectra.push(samples as Spectrum);
        }

        // For each registered dictionary
        for (index, dictionary) in dictionaries.iter_mut().enumerate() {
            // Each dictionary operates only in the fixed part of the spectrum
            // Selecting potentially interesting spectrum slice to check
            let frequencies = dictionary.get_bounds();
            let low  = (frequencies.0 / BASE_FREQUENCY).round() as usize;
            let high = (frequencies.1 / BASE_FREQUENCY).round() as usize;

            for fragment_index in 0 .. FRAGMENTS_PER_FRAME {
                let mut fragment_spectra = Vec::with_capacity(high - low);

                // We split the whole frame into several independent fragments each
                // having it's part of the time slice and own place in the dictionary
                let fragment_region = fragment_index*SLICES_PER_FRAGMENT .. (fragment_index + 1)*SLICES_PER_FRAGMENT;
                for spectrum in &spectra[fragment_region] {
                    // Spectrum slices for the particular fragment's frequency region
                    let slice: Spectrum = spectrum[low .. high + 1].iter().cloned().collect();
                    fragment_spectra.push(slice);
                }

                let fragment = Fragment::from_spectra(fragment_spectra);
                let fragment_key = dictionary.insert_fragment(fragment, 80);

                keys.push(fragment_key);
            }

            println!("frame {}, dictionary {}, fragments classified {}", frame_count, index, dictionary.len());
        }

        println!("processed frame {}\n", frame_count);
        frame_count += 1;
    }

    println!("\nCompleted.");

    (Glossary::from_dictionaries(detectors, dictionaries), keys)
}

pub fn reconstruct(filename: &str, glossary: &Glossary, keys: &KeyVec) {
    println!("Reconstructing {} from key vector of {} elements", filename, keys.len());

    let wav_header = hound::WavSpec {
        channels: 1,
        sample_rate: SAMPLE_RATE as u32,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };

    let mut writer = hound::WavWriter::create(filename, wav_header).unwrap();

    let plan = dft::Plan::new(dft::Operation::Backward, NUM_POINTS);
    let mut output = Vec::with_capacity(NUM_POINTS);
    let mut max_sample = 0.;

    let mut spectra = Vec::new();
    spectra.resize(SLICES_PER_FRAME, Spectrum::with_capacity(NUM_POINTS));

    let mut key_iter = keys.iter();

    for frame_index in 0 .. {
        println!("Reconstructing frame {}", frame_index);

        // Clearing from the previuos iteration
        output.clear();
        output.resize(NUM_POINTS, Cplx::default());

        for spectrum in spectra.iter_mut() {
            spectrum.clear();
            spectrum.resize(NUM_POINTS, Cplx::default());
        }

        // For each registered dictionary
        for (index, dictionary) in glossary.iter().enumerate() {
            // Each dictionary operates only in the fixed part of the spectrum
            // Selecting potentially interesting spectrum slice to check
            let frequencies = dictionary.get_bounds();
            let low  = (frequencies.0 / BASE_FREQUENCY).round() as usize;
            let high = (frequencies.1 / BASE_FREQUENCY).round() as usize;

            for fragment_index in 0 .. FRAGMENTS_PER_FRAME {
                let fragment_key = {
                    match key_iter.next() {
                        Some(option) => {
                            match option {
                                &Some(ref key) => key,
                                &None => continue
                            }
                        },

                        None => return
                    }
                };

                // Writing sub spectrum into it's place in the frame's spectra
                if let Some(fragment) = dictionary.find(fragment_key, 80) {
                    for (sub_index, sub_spectrum) in fragment.spectra().iter().enumerate() {
                        for (index, value) in sub_spectrum.iter().enumerate() {
                            spectra[fragment_index*SLICES_PER_FRAGMENT + sub_index][low + index] = *value;
                        }
                    }
                }
            }
        }

        // Performing IDFT on all spectrums of the frame
        for (slice_index, spectrum) in spectra.iter_mut().enumerate() {
            // Performing IFFT
            dft::transform(spectrum, &plan);

            // Applying Hanning window to the time domain values
            let normalized_length = 1. / (NUM_POINTS as f32 / 2.);
            for (sample_index, sample) in spectrum.iter_mut().take(NUM_POINTS / 2).enumerate() {
                let phase: f32 = Sample::from_sample(sample_index as f32 * normalized_length);
                let value: f32 = window::Hanning::at_phase(phase);

                *sample = *sample * (value / 2.);
            }

            // Mixing with output from other time slices
            let range = slice_index * SLICE_OFFSET .. slice_index * SLICE_OFFSET + NUM_POINTS / 2;
            for (sample_index, sample) in output[range].iter_mut().enumerate() {
                sample.re += spectrum[sample_index].re;
            }
        }

        // Normalizing output
        max_sample = max_sample.max(output.iter().max_by(|&x, &y| float_cmp(x.re, y.re, 0.00001)).unwrap().re);

        // Writing output to file
        for sample in &output {
            let amplitude = (i16::max_value() - 1000) as f32;
            writer.write_sample(((sample.re / max_sample) * amplitude * 0.75) as i16).unwrap();
        }

    }

    println!("done.");
}

pub fn dump_dictionary(filename: &str, dictionary: &Dictionary) {
    print!("Dumping dictionary {}... ", filename);

    let wav_header = hound::WavSpec {
        channels: 1,
        sample_rate: SAMPLE_RATE as u32,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };

    let mut writer = hound::WavWriter::create(filename, wav_header).unwrap();
    let plan = dft::Plan::new(dft::Operation::Backward, NUM_POINTS);

    let freqs = dictionary.get_bounds(); //FRAGMENT_WINDOW;
    let low  = (freqs.0 / BASE_FREQUENCY).round() as usize;
    let high = (freqs.1 / BASE_FREQUENCY).round() as usize;

    for (key, fragment) in dictionary.iter() {
//         println!("writing key {:?}\n", key);

        // Accumulated output mixed from all time slices
        let mut output = Vec::with_capacity(NUM_POINTS);
        output.resize(NUM_POINTS, Cplx::default());

        for (slice_index, samples) in fragment.spectra().iter()
            .map(|slice| {
                let mut spectrum = Vec::with_capacity(NUM_POINTS);

                // Preparing the full spectrum from it's fragment
                // by extending it with zeroes before and after
                spectrum.resize(low - 1, Cplx::default());
                spectrum.extend_from_slice(slice);
                spectrum.resize(NUM_POINTS, Cplx::default());

                // Performing IDFT
                dft::transform(&mut spectrum, &plan);

                // Applying Hanning window to the time domain values
                let normalized_length = 1. / (NUM_POINTS as f32 / 2.);
                for (sample_index, sample) in spectrum.iter_mut().take(NUM_POINTS / 2).enumerate() {
                    let phase: f32 = Sample::from_sample(sample_index as f32 * normalized_length);
                    let value: f32 = window::Hanning::at_phase(phase);

                    *sample = *sample * (value / 2.);
                }

                spectrum as Samples
            })
            .collect::<Vec<Samples>>()
            .iter()
            .enumerate()
        {
            // Mixing with output from other time slices
            let range = slice_index * SLICE_OFFSET .. slice_index * SLICE_OFFSET + NUM_POINTS / 2;
            for (sample_index, sample) in output[range].iter_mut().enumerate() {
                sample.re += samples[sample_index].re;
            }

            // Dumping individual spectrums
            {
                let max_sample = samples.iter().max_by(|&x, &y| float_cmp(x.re, y.re, 0.00001)).unwrap().re;
                for sample in &samples[.. NUM_POINTS / 2] {
                    let amplitude = (i16::max_value() - 1000) as f32;
                    writer.write_sample(((sample.re / max_sample) * amplitude * 0.75) as i16).unwrap();
                }

                // writing 100ms of silence between fragment slices
                for _ in 1 .. SAMPLE_RATE / 10 {
                    writer.write_sample(0).unwrap();
                }
            }
        }

        // Normalizing output
        let max_sample = output.iter().max_by(|&x, &y| float_cmp(x.re, y.re, 0.00001)).unwrap().re;

        // Writing output to file
        for sample in &output {
            let amplitude = (i16::max_value() - 1000) as f32;
            writer.write_sample(((sample.re / max_sample) * amplitude * 0.75) as i16).unwrap();
        }

        // Writing 1s silence between the fragments
        for _ in 1 .. SAMPLE_RATE {
            writer.write_sample(0).unwrap();
        }
    }

    println!("done.");
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
