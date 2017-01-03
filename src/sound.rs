
use num_complex::*;
use num_traits::Float;
// use std;

use dft;
use dft::{Operation, Plan};

pub use bit_vec::BitVec;

use hound;

// TODO Move to the Detector as individual field
const AMPLITUDE_DISPERSION_DB: f32 = 5.;

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

// pub fn filter_detectors(spectrum: &SpectrumSlice, detectors: &[Detector]) -> BitVec {
//     // This will hold resulting bit vector of the detectors activity mask
//     let mut result = BitVec::new();
//
//     filter_detectors_inplace(spectrum, detectors, &mut result);
//
//     result
// }

fn decibel(input: f32) -> f32 {
    20.0 * input.log10()
}

fn invert_decibel(input: f32) -> f32 {
    (10 as f32).powf(input / 20.0)
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

        println!("detector freq {}±{}, amp {}±{} dB, selected freq range is {} .. {}",
            detector.freq,
            detector.band,
            detector.amp,
            AMPLITUDE_DISPERSION_DB,
            freq(lo),
            freq(hi)
        );

        // Selecting the entry with the largest amplitude
        let (index, amplitude) = range
            .iter()
            .enumerate()
            .map(|(i, c)| (i, (c.norm() * 2.0) / NUM_POINTS as f32))
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

pub fn generate_file(filename: &str, detectors: &[Detector], mask: &BitVec) {
    let wav_header = hound::WavSpec {
        channels: 1,
        sample_rate: 44100,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };

    let mut writer = hound::WavWriter::create(filename, wav_header).unwrap();
    let mut mask_iter = mask.iter();

//     for frequency in (1 .. NUM_POINTS).map(|f| BASE_FREQUENCY * f as f32) {
//
//     }

    for chunk in 0 .. {
        let mut spectrum: Spectrum = Vec::with_capacity(NUM_POINTS);
        spectrum.resize(NUM_POINTS, Cplx::default());

        let mut num_taken = 0;

        for (freq, value) in mask_iter
            .clone()
            .take(detectors.len())
            .filter(|_| { num_taken += 1; true })
            .enumerate()
            .filter(|&(_, bit)| bit)
            .map(|(index, _)| {
                let frequency = detectors[index].freq;
                let amplitude = invert_decibel(detectors[index].amp) / 2.0;// * NUM_POINTS as f32;

                println!("chunk {}, active detector[{}]: freq {}, amplitude {} dB = {}",
                    chunk,
                    index,
                    frequency,
                    detectors[index].amp,
                    amplitude);

                // Phase was lost during analysis
                (frequency, Cplx::new(amplitude, 0.0))
            })
        {
            let bin = (freq / BASE_FREQUENCY).trunc() as usize;

            println!("writing bin[{}] = bin[{}] = {}", bin, NUM_POINTS-bin, value);

            // Spectrum has center symmetry
            spectrum[bin].re = value.re;
            spectrum[NUM_POINTS - bin].re = value.re;
        }

        if num_taken == 0 {
            break;
        }

        // Advance iterator to match cloned position
        mask_iter.nth(num_taken);

        // Synthesizing full range from the half by mirroring parts
        //let zero_iter = std::iter::once(Cplx::default());
//         let mut samples: Samples = spectrum
//             .iter()
//             .chain(vec![Cplx::default()].iter()) // FIXME Nyquist energy
//             .chain(spectrum.iter().rev())
//             .map(|x| *x)
//             .collect();

//         let mut samples = spectrum.clone();
//         samples.push(Cplx::default());
//         for &sample in spectrum.iter().rev() {
//             samples.push(sample);
//         }

        {
            let plan = Plan::new(Operation::Backward, NUM_POINTS);
            dft::transform(&mut spectrum, &plan);
        }

        let max_sample = spectrum.iter().max_by(|&x, &y| float_cmp(x.re, y.re, 0.00001)).unwrap().re;

        for sample in spectrum {
            let amplitude = (i16::max_value() - 100) as f32;
            writer.write_sample((sample.re / max_sample * amplitude) as i16).unwrap();
        }
    }
}
