
extern crate num_traits;
extern crate num_complex;
extern crate bit_vec;
extern crate byteorder;
extern crate dft;
extern crate clap;

use clap::{Arg, App, SubCommand};

use std::fs::File;
use std::cmp::Ordering;

use byteorder::{LittleEndian, ReadBytesExt/*, WriteBytesExt*/};
use dft::{Operation, Plan/*, c64*/};

mod sound;

use sound::*;

fn main() {
    let options = App::new("Semantic sound processor")
        .version("0.1")
        .author("Dmitriy Kashitsyn <korvin@deeptown.org>")
        .about("Semantic sound processor based on Alex Redozubov's pattern wave theory")
        .arg(Arg::with_name("input")
            .short("i")
            .long("input")
            .help("Sets the input file to use")
            .required(true)
            .takes_value(true))
        .get_matches();

    let input_filename = options.value_of("input").unwrap();
    let mut file = File::open(input_filename).unwrap();

    let mut data = Vec::new();
    while let Ok(input) = file.read_i16::<LittleEndian>() {
        data.push(sound::Cplx::new(input as f32, 0.0));
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

//         for freq in 1 .. 20 {
//             detectors.push(Detector::new(freq as f32 * 100.0, 50.0, 8000.0));
//             detectors.push(Detector::new(freq as f32 * 100.0, 50.0, 12000.0));
//             detectors.push(Detector::new(freq as f32 * 100.0, 50.0, 16000.0));
//             detectors.push(Detector::new(freq as f32 * 100.0 , 50.0, 20000.0));
//         }

        // Populating detectors from 0Hz to ~1KHz with 100Hz selectivity (±50 Hz)
        for i in 1 .. 12 {
            let freq = BASE_FREQUENCY * 2.0 * i as f32;

            detectors.push(Detector::new(freq, 50.0, 8000.0));
            detectors.push(Detector::new(freq, 50.0, 16000.0));
            detectors.push(Detector::new(freq, 50.0, 12000.0));
            detectors.push(Detector::new(freq, 50.0, 20000.0));
        }

        // Populating detectors from ~1Hz to 3KHz with 500Hz selectivity
        for i in 0 .. 10 {
//             let freq = 990.0 + BASE_FREQUENCY * i as f32;
//             let band = 24.7 * (4.37 * freq as f32 + 1.0);
            let freq = 990.0 + 4.0 * BASE_FREQUENCY * i as f32;

            detectors.push(Detector::new(freq, 500.0, 8000.0));
            detectors.push(Detector::new(freq, 500.0, 16000.0));
            detectors.push(Detector::new(freq, 500.0, 12000.0));
            detectors.push(Detector::new(freq, 500.0, 20000.0));
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
