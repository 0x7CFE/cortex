
extern crate num_traits;
extern crate num_complex;
extern crate bit_vec;
// extern crate byteorder;
extern crate dft;
extern crate clap;
extern crate hound;

use clap::{Arg, App, SubCommand};

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
    let mut detectors = Vec::new();

//     detectors.push(Detector::new(990.0, 50.0, -5.0));
//     detectors.push(Detector::new(990.0, 50.0, -15.0));
//     detectors.push(Detector::new(990.0, 50.0, -25.0));
//     detectors.push(Detector::new(990.0, 50.0, -35.0));

    // Populating detectors from 0Hz to ~1KHz with 100Hz selectivity (±50 Hz)
    for i in 1 .. 12 {
        let freq = BASE_FREQUENCY * 2.0 * i as f32;

        detectors.push(Detector::new(freq, 50.0, -5.0));
        detectors.push(Detector::new(freq, 50.0, -15.0));
        detectors.push(Detector::new(freq, 50.0, -25.0));
        detectors.push(Detector::new(freq, 50.0, -35.0));
    }

    // Populating detectors from ~1Hz to 3KHz with 500Hz selectivity
    for i in 0 .. 10 {
        let freq = 990.0 + 4.0 * BASE_FREQUENCY * i as f32;

        detectors.push(Detector::new(freq, 500.0, -5.0));
        detectors.push(Detector::new(freq, 500.0, -15.0));
        detectors.push(Detector::new(freq, 500.0, -25.0));
        detectors.push(Detector::new(freq, 500.0, -35.0));
    }

    let mask = analyze_file(input_filename, &detectors);
    //let mask = sound::filter_detectors(result, &detectors);
    println!("{} detectors, mask size {} : {:?}", detectors.len(), mask.len(), mask);

    /*let plan = Plan::new(Operation::Backward, NUM_POINTS);
    dft::transform(&mut data[..NUM_POINTS], &plan);

    let mut file = File::create("output.pcm").unwrap();

    for c in data.iter()/*.take(10)*/ {
        let n = c.re / NUM_POINTS as f64;
        let v = n.round() as i16;

//         println!("writing c {} -> {} -> {}", c, n, v);
        file.write_i16::<LittleEndian>(v).unwrap();
    }*/

}
