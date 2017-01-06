
use num_complex::*;
use num_traits::Float;
use std::f32::consts::PI;

use dft;
use dft::{Operation, Plan};

pub use bit_vec::BitVec;

use hound;

// TODO Move to the Detector as individual field
const AMPLITUDE_DEVIATION_DB: f32 = 5.;
const PHASE_DEVIATION_DB: f32 = PI / 4.;

pub type Cplx = Complex<f32>;

use std::cmp::Ordering;
pub type Samples  = Vec<Cplx>;
pub type Spectrum = Vec<Cplx>;
pub type SpectrumSlice = [Cplx];

pub const SAMPLE_RATE:     usize = 96000; // 44100;
pub const NUM_POINTS:      usize = 8192;
pub const BASE_FREQUENCY:  f32 = (SAMPLE_RATE as f32) / (NUM_POINTS as f32);

#[derive(Copy, Clone, Debug, Default)]
pub struct Detector {
    freq:  f32, // base detector frequency
    band:  f32, // frequency range
    amp:   f32, // amplitude
    phase: f32, // phase
}

impl Detector {
    pub fn new(freq: f32, band: f32, amp: f32, phase: f32) -> Detector {
        Detector {
            freq: freq,
            band: band,
            amp: amp,
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

fn decibel(input: f32) -> f32 {
    20.0 * input.log10()
}

fn invert_decibel(input: f32) -> f32 {
    (10 as f32).powf(input / 20.0)
}

pub fn filter_detectors_inplace(spectrum: &SpectrumSlice, detectors: &[Detector], result: &mut BitVec) {
    // Iterating through all detectors filtering out activity
    for detector in detectors {
        // Each detector operates only in the fixed part of the spectrum
        // Selecting potentially interesting spectrum slice to check

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
        let amp_match   = (decibel(amplitude).abs() - detector.amp.abs()).abs() < AMPLITUDE_DEVIATION_DB;
        let phase_match = ((phase + PI) - (detector.phase + PI)).abs() < PHASE_DEVIATION_DB;
        let is_active   = amp_match && phase_match;

        println!("signal amplitude {} → {} dB, phase {}{}\n",
            amplitude,
            decibel(amplitude),
            phase,
            if is_active { ", **MATCH**" } else { "" }
        );

        result.push(is_active);
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
            //.max_by(|&(_, x, _), &(_, y, _)| x.partial_cmp(&y).unwrap_or(Ordering::Less))
            .max_by(|&(_, x, _), &(_, y, _)| float_cmp(x, y, 0.00001))
            .unwrap();

        // Treating detector as active if max amplitude lays within detector's selectivity range
        let amp_match   = (decibel(amplitude).abs() - detector.amp.abs()).abs() < AMPLITUDE_DEVIATION_DB;
        let phase_match = ((phase + PI) - (detector.phase + PI)).abs() < PHASE_DEVIATION_DB;
        let is_active   = amp_match && phase_match;

        println!("signal frequency {}, amplitude {} → {} dB, phase {}{}\n",
            freq(lo+index),
            amplitude,
            decibel(amplitude),
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
            .filter(|_| { total_samples += 1; true })
            .take(NUM_POINTS)
            .map(|s| Cplx::new(s.unwrap() as f32 / i16::max_value() as f32, 0.0))
            .collect();

        if samples.len() < NUM_POINTS / 2 {
            println!("*** End of input, samples processed: {}, chunks {}", total_samples, chunk - 1);
            break;
        }

        println!("*** Chunk {}", chunk);
        samples.resize(NUM_POINTS, Cplx::default());

        let plan = Plan::new(Operation::Forward, NUM_POINTS);
        dft::transform(&mut samples, &plan);

        filter_detectors_inplace(&samples[..NUM_POINTS/2 - 1], detectors, &mut result);
        println!("total samples so far {}", total_samples);
    }

    result
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
                let amplitude = invert_decibel(detectors[index].amp) / 2.0;// * NUM_POINTS as f32;
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

        let plan = Plan::new(Operation::Backward, NUM_POINTS);
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
