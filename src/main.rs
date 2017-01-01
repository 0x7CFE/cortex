#![feature(iter_max_by)]

extern crate num_traits;
extern crate num_complex;
extern crate bit_vec;
extern crate byteorder;
extern crate dft;

use std::fs::File;
use std::cmp::Ordering;

use byteorder::{LittleEndian, ReadBytesExt/*, WriteBytesExt*/};
use dft::{Operation, Plan/*, c64*/};

mod sound;

use sound::*;

fn main() {
    let mut file = File::open("input.pcm").unwrap();
    let mut data = Vec::new();

    while let Ok(input) = file.read_i16::<LittleEndian>() {
        data.push(sound::Cplx::new(input as f64, 0.0));
    }

    // Filling the remaining part of the window with zeroes
    if data.len() < NUM_POINTS {
        data.resize(NUM_POINTS, sound::Cplx::new(0.0, 0.0));
    }

    let plan = Plan::new(Operation::Forward, NUM_POINTS);
    dft::transform(&mut data[..NUM_POINTS], &plan);

    {
        // Fourier transform of real input is symmetrical, we need only the half
        let result = &data[..NUM_POINTS/2-1];

        let (index, magnitude) = result.iter().enumerate().map(|(i, c)| (i, c.norm())).max_by(|&(_, x), &(_, y)| {
            if (x - y).abs() < 0.00001 { return Ordering::Equal; }

            if x < y {
                Ordering::Less
            } else {
                Ordering::Greater
            }
        }).unwrap();

//         let magnitude = magnitude / NUM_POINTS as f64;

        let frequency = (index as f64) * (SAMPLE_RATE as f64) / (NUM_POINTS as f64);
        println!("global max frequency is {}, max magnitude is {} and it's index is {}", frequency, magnitude, index);

        let mut detectors = Vec::new();

        // Populating detectors from 100Hz to 2KHz with 100Hz selectivity
        for freq in 1 .. 20 {
            detectors.push(Detector::new(freq * 100, 50, 8000));
            detectors.push(Detector::new(freq * 100, 50, 12000));
            detectors.push(Detector::new(freq * 100, 50, 16000));
            detectors.push(Detector::new(freq * 100, 50, 20000));

//             detectors.push(Detector::new(freq * 150, 50, 0000));
//             detectors.push(Detector::new(freq * 150, 50, 4000));
//             detectors.push(Detector::new(freq * 150, 50, 8000));
//             detectors.push(Detector::new(freq * 150, 50, 16000));
        }
/*
        let detectors = vec![
            Detector::new(387, 50, 19000),
            Detector::new(1000, 50, 19000),
            Detector::new(2000, 50, 19000),
        ];*/

        let mask = sound::filter_detectors(result, &detectors);
        println!("detector activity mask is {:?}", mask);

    }

    /*let plan = Plan::new(Operation::Backward, NUM_POINTS);
    dft::transform(&mut data[..NUM_POINTS], &plan);

    let mut file = File::create("output.pcm").unwrap();

    for c in data.iter()/*.take(10)*/ {
        let n = c.re / NUM_POINTS as f64;
        let v = n.round() as i16;

//         println!("writing c {} -> {} -> {}", c, n, v);
        file.write_i16::<LittleEndian>(v).unwrap();
    }*/


//     calculate the magnitude of each DFT output bin: magnitude = sqrt(re*re+im*im)
//     find the bin with the largest magnitude, call its index i_max.
//     calculate the equivalent frequency of this bin: freq = i_max * Fs / N, here Fs = sample rate (Hz) and N = no of points in FFT.

}

// Расположение и структура формант (соотношение их амплитуд, добротностей и др.) являются основным признаком
// различения гласных звуков и в значительной степени согласных. При обработке звука во внутреннем ухе происходит
//  спектральный анализ сигнала, при этом разрешающая способность этого анализа зависит от ширины «слуховых»
//  фильтров на базилярной мембране («критических полос слуха»). Ширина этих фильтров зависит от частоты:
//  примерно до 1 кГц она равна 100 Гц, выше – меняется по закону:
// Δ=24,7(4,37 fср + 1), где fср – центральная частота в каждой полосе.
//
