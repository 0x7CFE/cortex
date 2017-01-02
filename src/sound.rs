

use num_complex::Complex;
use num_traits::{Num, Float};

use dft;
use dft::{Operation, Plan};

pub use bit_vec::BitVec;

use hound;

// TODO Move to the Detector as individual field
const AMPLITUDE_DISPERSION_DB: f32 = 10.0;

pub type Cplx = Complex<f32>;

use std::cmp::Ordering;
pub type Samples  = Vec<Cplx>;
pub type Spectrum = Vec<Cplx>;
pub type SpectrumSlice = [Cplx];

pub const SAMPLE_RATE:     usize = 44100;
pub const NUM_POINTS:      usize = 1024;
pub const BASE_FREQUENCY:  f32 = SAMPLE_RATE as f32 / NUM_POINTS as f32;

#[derive(Copy, Clone, Debug, Default)]
pub struct Detector {
    freq: f32, // base detector frequency
    band: f32, // frequency range
    amp:  f32, // amplitude
//     phase: i16,
}

impl Detector {
    pub fn new(freq: f32, band: f32, amp: f32/*, phase: i16*/) -> Detector {
        Detector {
            freq: freq,
            band: band,
            amp: amp,

//             phase: phase
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

pub fn filter_detectors(spectrum: &SpectrumSlice, detectors: &[Detector]) -> BitVec {
    // This will hold resulting bit vector of the detectors activity mask
    let mut result = BitVec::new();

    filter_detectors_inplace(spectrum, detectors, &mut result);

    result
}

fn decibel(input: f32) -> f32 {
    20.0 * input.log10()
}

pub fn filter_detectors_inplace(spectrum: &SpectrumSlice, detectors: &[Detector], result: &mut BitVec) {
//     println!("base frequency is {}, amp response {}, spectrum len {}\n",
//         BASE_FREQUENCY,
//         AMPLITUDE_DISPERSION,
//         spectrum.len());

    // Iterating through all detectors filtering out activity
    for detector in detectors {
        // Each detector operates only in the fixed part of the spectrum
        // Selecting potentially interesting spectrum slice to check
        let lo = ((detector.freq - detector.band) / BASE_FREQUENCY).round() as usize;
        let hi = ((detector.freq + detector.band) / BASE_FREQUENCY).round() as usize;

        if lo > spectrum.len() - 1 || hi > spectrum.len() - 1 {
            println!("invalid detector freq {}, band {}", detector.freq, detector.band);
            break;
        }

        let range =
            if hi + 1 < spectrum.len() {
                &spectrum[lo .. hi + 1]
            } else {
                &spectrum[lo .. spectrum.len() - 1]
            };

        println!("detector freq {}±{}, amp {}, selected freq range is {} .. {}",
            detector.freq,
            detector.band,
            detector.amp,
            freq(lo),
            freq(hi)
        );

        // Selecting the entry with the largest amplitude
        let (index, amplitude) = range
            .iter()
            .enumerate()
            .map(|(i, c)| (i, c.norm() * 2.0 / NUM_POINTS as f32))
            .max_by(|&(_, x), &(_, y)| float_cmp(x, y, 0.00001))
            .unwrap();

        // Treating detector as active if max amplitude lays within detector's selectivity range
        let is_active = (decibel(amplitude).abs() - detector.amp.abs()).abs() < AMPLITUDE_DISPERSION_DB;

        println!("signal frequency {}, amplitude {} → {} dB {}\n",
            freq(lo+index),
            amplitude,
            decibel(amplitude),
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

    // Reading file by chunks of NUM_POINTS, then processing
    for i in 0 .. {
        let mut samples: Samples = reader
            .samples::<i16>()
            .take(NUM_POINTS)
            .map(|s| Cplx::new(s.unwrap() as f32 / i16::max_value() as f32, 0.0))
            .collect();

        if samples.len() < NUM_POINTS {
            break;
        }

        println!("*** Chunk {}", i);

        let plan = Plan::new(Operation::Forward, NUM_POINTS);
        dft::transform(&mut samples, &plan);

        filter_detectors_inplace(&samples[..NUM_POINTS/2 - 1], detectors, &mut result);
    }

    result
}
