

use num_complex::Complex;
use num_traits::Float; //*, FloatConst, One*/};
pub use bit_vec::BitVec;

const AMPLITUDE_DISPERSION: f32 = 2000.0 * 1.0;

pub type Cplx = Complex<f32>;

use std::cmp::Ordering;
pub type Spectrum = Vec<Complex<f32>>;
pub type SpectrumSlice = [Complex<f32>];

pub const SAMPLE_RATE:     usize = 44100;
pub const NUM_POINTS:      usize = 1024;
pub const BASE_FREQUENCY:  f32 = SAMPLE_RATE as f32 / NUM_POINTS as f32;

#[derive(Copy, Clone, Debug, Default)]
pub struct Detector {
    freq: f32,
    band: f32,
    amp:  f32,
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

    /*pub fn match_freq(&self, freq: i16) -> bool {
        (self.freq - freq).abs() < self.band
    }

    pub fn match_amp(&self, amp: i16) -> bool {
        (self.amp - amp).abs() < 4096
    }*/
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

    println!("base frequency is {}, amp response {}, spectrum len {}\n",
        BASE_FREQUENCY,
        AMPLITUDE_DISPERSION,
        spectrum.len());

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
            .map(|(i, c)| (i, c.norm()))
            .max_by(|&(_, x), &(_, y)| float_cmp(x, y, 0.00001))
            .unwrap();

        // Treating detector as active if max amplitude lays within detector's selectivity range
        let detector_amplitude = detector.amp as f32;
        let is_active = (amplitude.abs() - detector_amplitude).abs() < AMPLITUDE_DISPERSION;

        println!("signal frequency {}, amplitude {}{}\n",
            freq(lo+index),
            amplitude,
            if is_active { ", **MATCH**" } else { "" }
        );

        result.push(is_active);
    }

    result
}

