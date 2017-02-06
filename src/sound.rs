
use std::sync::Arc;
use std::f32::consts::PI;
use std::io::{Read};
use dft;

use hound;
use sample::*;
use sample::window::Type;

use num_complex;
use num_traits::Float;

use itertools::Itertools;
use iter::{CollectChunks, CollectChunksExt};

use rayon::prelude::*;

use memory::*;

// TODO Move to the Detector as individual field
pub const AMPLITUDE_DEVIATION_DB: f32 = 5.;

pub type Cplx = ::num_complex::Complex<f32>;

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

#[derive(Copy, Clone, Debug, Default, Serialize, Deserialize)]
pub struct Detector {
    pub freq:  f32, // base detector frequency
    pub band:  f32, // frequency range
    pub amp:   f32, // amplitude in dB
    pub phase: f32, // phase
    pub phase_range: f32,
}

pub fn detector_freq(index: usize) -> f32 {
    let ideal_freq = 15. + 5. * index as f32 + ((index as f32 - 5.) / 16.).exp();
    let fft_freq = (ideal_freq / BASE_FREQUENCY).trunc() * BASE_FREQUENCY;

    fft_freq
}

impl Detector {
    pub fn new(freq: f32, band: f32, amp: f32, phase: f32, phase_range: f32) -> Detector {
        Detector {
            freq:  freq,
            band:  band,
            amp:   amp,
            phase: phase,
            phase_range: phase_range,
        }
    }
}

fn freq(index: usize) -> f32 {
    (index as f32) * (SAMPLE_RATE as f32) / (NUM_POINTS as f32)
}

// TODO Rewrite using https://crates.io/crates/float-cmp
pub fn float_cmp<F: Float>(x: F, y: F, epsilon: F) -> Ordering {
    if (x - y).abs() < epsilon { return Ordering::Equal; }

    if x < y {
        Ordering::Less
    } else {
        Ordering::Greater
    }
}

const SLICES_PER_FRAME:    usize = 16;
const FRAGMENTS_PER_FRAME: usize = 4;

const SLICE_OFFSET:        usize = (NUM_POINTS / 2) / SLICES_PER_FRAME;
const SLICES_PER_FRAGMENT: usize = SLICES_PER_FRAME / FRAGMENTS_PER_FRAME;

pub type KeyVec = Vec<Option<FragmentKey>>;

pub fn build_detectors() -> Vec<Detector> {
    let mut detectors = Vec::new();

    for i in 10 .. 141 {
        let freq = detector_freq(i);
        let band = 2. * (detector_freq(i+1) - freq);

        let phase_count = match i {
            _ if detector_freq(i) < 100.  => 30,
            _ if detector_freq(i) < 300.  => 25,
            _ if detector_freq(i) < 500.  => 20,
            _ if detector_freq(i) < 1000. => 10,
            _ if detector_freq(i) < 2000. => 8,
            _ => 4,
        };

        let phase_range = 6. * PI / phase_count as f32;

        for phase in (0 .. phase_count).map(|n| 2. * PI * n as f32 / phase_count as f32 - PI)
        {
            detectors.push(Detector::new(freq, band, -5.,  phase, phase_range));
            detectors.push(Detector::new(freq, band, -15., phase, phase_range));
            detectors.push(Detector::new(freq, band, -25., phase, phase_range));
            detectors.push(Detector::new(freq, band, -35., phase, phase_range));
        }
    }

    detectors
}

pub fn build_bounds() -> Vec<(f32, f32)> {
    (1 .. 14).into_iter()
        .map(|i| (detector_freq(i*10), detector_freq((i+1)*10)) )
        .collect()
}

pub fn build_dictionaries(bounds: &[(f32, f32)]) -> Vec<Dictionary> {
    bounds.iter()
        .map(|&(low, high)| Dictionary::new(low, high))
        .collect()
}

fn build_spectra(first_chunk: &Samples, second_chunk: &Samples, plan: &dft::Plan<f32>) -> Vec<Spectrum> {
    (0 .. SLICES_PER_FRAME).into_par_iter()
        .map( |slice| {
            // Collecting samples
            let first_range  = slice * SLICE_OFFSET ..;
            let second_range = .. slice * SLICE_OFFSET;

            let mut samples: Samples = first_chunk[first_range].iter()
                .chain(second_chunk[second_range].iter())
                .cloned()
                .collect();

            // Zero padding to the full frame length
            samples.resize(NUM_POINTS, Cplx::default());

            // Performing FFT
            dft::transform(&mut samples, plan);

            samples
        })
        .collect()
}

pub fn analyze_file(filename: &str) -> KeyVec {
    println!("Generating keys for {}...", filename);

    // This will hold resulting bit vector of the detectors activity mask
    let mut result = KeyVec::new();

    let detectors = build_detectors();
    let bounds = build_bounds();

    // Opening wave file for reading
    let mut reader = hound::WavReader::open(filename).unwrap();
    let plan = dft::Plan::new(dft::Operation::Forward, NUM_POINTS);

    let mut total_samples = 0;

    for (first, second) in reader
        .samples::<i16>()
        .map(|s| Cplx::new(s.unwrap() as f32 / i16::max_value() as f32, 0.0))
        .collect_chunks::<Samples>(NUM_POINTS / 2)
        .map(Arc::new) // Wrap each chunk into Arc for cheap cloning
        .tuple_windows()
    {
        if first.len() < NUM_POINTS / 2 || second.len() < NUM_POINTS / 2 {
            break;
        }

        // Collecting overlapping spectrums spanning the whole frame
        let spectra: Vec<_> = build_spectra(&first, &second, &plan);

        // For each registered dictionary
        for (index, bound) in bounds.iter().enumerate() {
            // Each dictionary operates only in the fixed part of the spectrum
            // Selecting potentially interesting spectrum slice to check
            let low  = (bound.0 / BASE_FREQUENCY).round() as usize;
            let high = (bound.1 / BASE_FREQUENCY).round() as usize;

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
                let key = fragment.get_key(*bound, &detectors);

                // Appending fragment key to result
                if key.bits_set() != 0 {
                    result.push(Some(key))
                } else {
                    result.push(None)
                }
            }
        }
    }

    println!("\nCompleted.");

    result
}

pub fn build_glossary(filename: &str, similarity: usize) -> Glossary {
    println!("Building glossary from {} ... ", filename);

    let detectors = build_detectors();
    let mut dictionaries = build_dictionaries(&build_bounds());

    let plan = dft::Plan::new(dft::Operation::Forward, NUM_POINTS);
    let mut reader = hound::WavReader::open(filename).unwrap();

    let mut frame_count = 0;

    for (first, second) in reader
        .samples::<i16>()
        .map(|s| Cplx::new(s.unwrap() as f32 / i16::max_value() as f32, 0.0))
        .collect_chunks::<Samples>(NUM_POINTS / 2)
        .map(Arc::new) // Wrap each chunk into Arc for cheap cloning
        .tuple_windows()
    {
        if first.len() < NUM_POINTS / 2 || second.len() < NUM_POINTS / 2 {
            break;
        }

        // Collecting overlapping spectrums spanning the whole frame
        let spectra: Vec<_> = build_spectra(&first, &second, &plan);

        // For each registered dictionary
        dictionaries.par_iter_mut().enumerate().for_each( |(index, dictionary)| {
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
                let key = dictionary.insert_fragment(fragment, &detectors, similarity);
            }

            println!("frame {}, dictionary {}, fragments classified {}", frame_count, index, dictionary.len());
        });

        println!();
        frame_count += 1;
    }

    println!("\nCompleted.");

    Glossary::new(detectors, dictionaries)
}

pub fn reconstruct(filename: &str, glossary: &Glossary, keys: &KeyVec, similarity: usize) {
    println!("Reconstructing {} from key vector of {} elements", filename, keys.len());

    let wav_header = hound::WavSpec {
        channels: 1,
        sample_rate: SAMPLE_RATE as u32,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };

    let mut writer = hound::WavWriter::create(filename, wav_header).unwrap();

    let plan = dft::Plan::new(dft::Operation::Backward, NUM_POINTS);
    let mut max_sample = 0.6;

    let mut output = Vec::with_capacity(NUM_POINTS);
    output.resize(NUM_POINTS, Cplx::default());

    let mut spectra = Vec::new();
    spectra.resize(SLICES_PER_FRAME, Spectrum::with_capacity(NUM_POINTS));

    let mut key_iter = keys.iter();

    for frame_index in 0 .. {
        println!("Reconstructing frame {}", frame_index);

        // Left shift of the output by the half
        output = output.split_off(NUM_POINTS / 2);
        output.resize(NUM_POINTS, Cplx::default());

        // Clearing from the previuos iteration
        for spectrum in &mut spectra {
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
                            match *option {
                                Some(ref key) => key,
                                None => continue
                            }
                        },

                        None => return
                    }
                };

                // Writing sub spectrum into it's place in the frame's spectra
                if let Some(fragment) = dictionary.find(fragment_key, similarity) {
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
        for sample in &output[ .. NUM_POINTS / 2] {
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
