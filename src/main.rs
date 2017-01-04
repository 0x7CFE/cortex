
extern crate num_traits;
extern crate num_complex;
extern crate bit_vec;
extern crate dft;
extern crate clap;
extern crate hound;

use clap::{Arg, App};

mod sound;
use sound::Detector;

fn detector_freq(index: usize) -> f32 {
    15. + 5. * index as f32 + ((index as f32 - 5.) / 16.).exp()
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
        let band = 2. * (freq - detector_freq(i-1));

        detectors.push(Detector::new(freq, band, -5.));
        detectors.push(Detector::new(freq, band, -15.));
        detectors.push(Detector::new(freq, band, -25.));
        detectors.push(Detector::new(freq, band, -35.));

        println!("detector[{:3}]\tfreq {:.2},\tband {:.2}", i, freq, band);
    }

//     detectors.push(Detector::new(516.7969, 50.0, -12.0));
//     detectors.push(Detector::new(990.52734, 50.0, -5.0));
//     detectors.push(Detector::new(1050.52734, 50.0, -5.0));

//    detectors.push(Detector::new(1990.0, 50.0, -15.0));
//     detectors.push(Detector::new(990.0, 50.0, -25.0));
//     detectors.push(Detector::new(990.0, 50.0, -35.0));

    let mask = sound::analyze_file(input_filename, &detectors);
    println!("{} detectors, mask size {}", detectors.len(), mask.len());

    for (i, bit) in mask.iter().enumerate() {
        if i % detectors.len() == 0 {
            println!();
            print!("{:04} : ", i / detectors.len() + 1);
        }

        if i % (detectors.len() / 4) == 0 {
            print!("\n        ");
        }

        print!("{}", bit as i8);
    }

    print!("\n\n");

    if let Some(output_filename) = options.value_of("output") {
        sound::generate_file(output_filename, &detectors, &mask);
    }
}
