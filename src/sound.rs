

use num_complex::Complex;
use num_traits::Float; //*, FloatConst, One*/};
pub use bit_vec::BitVec;

const AMPLITUDE_DISPERSION: f64 = 2000.0 * 1.0;

pub type Cplx = Complex<f64>;

use std::cmp::Ordering;
pub type Spectrum = Vec<Complex<f64>>;
pub type SpectrumSlice = [Complex<f64>];

pub const SAMPLE_RATE:     usize = 44100;
pub const NUM_POINTS:      usize = 1024;
pub const BASE_FREQUENCY:  usize = SAMPLE_RATE / NUM_POINTS;

#[derive(Copy, Clone, Hash, Debug, Default)]
pub struct Detector {
    frequency: i16,
    dispersion: i16,
    amplitude: i16,
//     phase: i16,
}

impl Detector {
    pub fn new(freq: i16, disp: i16, amp: i16/*, phase: i16*/) -> Detector {
        Detector {
            frequency: freq,
            dispersion: disp,
            amplitude: amp,
//             phase: phase
        }
    }

    /*pub fn match_freq(&self, freq: i16) -> bool {
        (self.frequency - freq).abs() < self.dispersion
    }

    pub fn match_amp(&self, amp: i16) -> bool {
        (self.amplitude - amp).abs() < 4096
    }*/
}

fn freq(index: usize) -> f64 {
    (index as f64) * (SAMPLE_RATE as f64) / (NUM_POINTS as f64)
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

    println!("base frequency is {}, amplitude response {}, spectrum len {}\n",
        BASE_FREQUENCY,
        AMPLITUDE_DISPERSION,
        spectrum.len());

    // Iterating through all detectors filtering out activity
    for detector in detectors {
        // Each detector operates only in the fixed part of the spectrum
        // Selecting potentially interesting spectrum slice to check
        let lo = ((detector.frequency - detector.dispersion) / BASE_FREQUENCY as i16) as usize;
        let hi = ((detector.frequency + detector.dispersion) / BASE_FREQUENCY as i16) as usize;
        let range = &spectrum[lo .. hi + 1];

        println!("detector freq {}±{}, amp {}, selected freq range is {} .. {}",
            detector.frequency,
            detector.dispersion,
            detector.amplitude,
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
        let detector_amplitude = detector.amplitude as f64;
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

