
#![feature(btree_range, collections_bound)]
#![feature(conservative_impl_trait)]

extern crate num_traits;
extern crate num_complex;
extern crate bit_vec;
extern crate dft;
extern crate clap;
extern crate hound;

#[macro_use] extern crate itertools;

use clap::{Arg, App};

mod sound;
mod memory;
mod iter;

use sound::Detector;

use std::f32::consts::PI;

fn detector_freq(index: usize) -> f32 {
    let ideal_freq = 15. + 5. * index as f32 + ((index as f32 - 5.) / 16.).exp();
    let fft_freq = (ideal_freq / sound::BASE_FREQUENCY).trunc() * sound::BASE_FREQUENCY;

    fft_freq
}

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
        .arg(Arg::with_name("output")
            .short("o")
            .long("output")
            .help("Sets the output file to rewrite")
            .required(false)
            .takes_value(true))
        .get_matches();

    let input_filename = options.value_of("input").unwrap();
    let mut detectors = Vec::new();

    // Populating detectors from 0Hz to ~1KHz with ~2683Hz
    for i in 1 .. 129 {
        let freq = detector_freq(i);
        let band = 2. * (detector_freq(i+1) - freq);

//         for phase in vec![ -PI, -3.*PI/4., -PI/2., -PI/4., 0., PI/4., PI/2., 3.*PI/4., PI]
        for phase in vec![ -PI, -PI/2., 0., PI/2., PI ]
        {
            detectors.push(Detector::new(freq, band, -5.,  phase));
            detectors.push(Detector::new(freq, band, -15., phase));
            detectors.push(Detector::new(freq, band, -25., phase));
            detectors.push(Detector::new(freq, band, -35., phase));
            //detectors.push(Detector::new(freq, band, -45., phase));
        }

        println!("detector[{:3}]\tfreq {:.2},\tband {:.2}", i, freq, band);
    }

    let dictionary = sound::build_dictionary(input_filename, &detectors);

    for (key, value) in dictionary.iter() {
        //println!("{:?} -> {:?}\n", key, value);
        println!("{:?}\n", key);
    }

    println!("{} fragments total", dictionary.len());

    if let Some(output_filename) = options.value_of("output") {
        sound::dump_dictionary(output_filename, &dictionary);
    }


    /*let mask = sound::analyze_file(input_filename, &detectors);
    println!("{} detectors, mask size {}", detectors.len(), mask.len());

    for (i, bit) in mask.iter().enumerate() {
        if i % detectors.len() == 0 {
            println!();
            print!("{:04} : ", i / detectors.len() + 1);
        }

        if i % (detectors.len() / 16) == 0 {
            print!("\n        ");
        }

        print!("{}", if bit { "!" } else { "." } );
    }

    print!("\n\n");

    if let Some(output_filename) = options.value_of("output") {
        sound::generate_file(output_filename, &detectors, &mask);
    }*/
}
