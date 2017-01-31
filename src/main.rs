
#![feature(btree_range, collections_bound)]
#![feature(conservative_impl_trait)]

extern crate num_traits;
extern crate num_complex;
extern crate bit_vec;
extern crate dft;
extern crate clap;
extern crate hound;
extern crate sample;
extern crate rayon;

#[macro_use] extern crate itertools;
#[macro_use] extern crate serde_derive;
extern crate serde;

// extern crate serde_json;
// extern crate bincode;
// use bincode::SizeLimit;
// use bincode::serde::{serialize_into, deserialize_from};

extern crate serde_cbor;
use serde_cbor::ser::to_writer;
use serde_cbor::de::from_reader;

use clap::{Arg, App};

mod sound;
mod memory;
mod iter;

use sound::Detector;

use sound::*;
use memory::*;

use std::f32::consts::PI;
use std::fs::File;

fn main() {
    let options = App::new("Semantic sound processor")
        .version("0.1")
        .author("Dmitriy Kashitsyn <korvin@deeptown.org>")
        .about("Semantic sound processor based on Alex Redozubov's pattern wave theory")
        .arg(Arg::with_name("input")
            .short("i")
            .long("input")
            .help("Sets the input file to use")
            .required(false)
            .takes_value(true))
//         .arg(Arg::with_name("output")
//             .short("o")
//             .long("output")
//             .help("Sets the output file to rewrite")
//             .required(false)
//             .takes_value(true))
        .arg(Arg::with_name("similarity")
            .long("similarity")
            .help("Similarity coefficient in percents")
            .required(false)
            .takes_value(true))
        .arg(Arg::with_name("build-glossary")
            .long("build-glossary")
            .help("Sets the glossary file name to build")
            .required(false)
            .takes_value(true))
        .arg(Arg::with_name("use-glossary")
            .long("use-glossary")
            .help("Sets the glossary file to use")
            .required(false)
            .takes_value(true))
        .arg(Arg::with_name("use-key")
            .long("use-key")
            .help("Sets the key file to use")
            .required(false)
            .takes_value(true))
        .arg(Arg::with_name("write-key")
            .long("write-key")
            .help("Sets the key file to write")
            .required(false)
            .takes_value(true))
        .arg(Arg::with_name("reconstruct")
            .long("reconstruct")
            .help("Sets the reconstruction file name")
            .required(false)
            .takes_value(true))
        .get_matches();

    let similarity =
        if let Some(value) = options.value_of("similarity") {
            value.parse().unwrap()
        } else {
            50
        };

    let (glossary, keys) = {
        if let Some(glossary_filename) = options.value_of("build-glossary") {
            let input_filename = options.value_of("input").unwrap();
            let (glossary, keys) = sound::build_glossary(input_filename, similarity);

            println!("Writing glossary file");
            let mut glossary_file = File::create(glossary_filename).unwrap();
//            serialize_into(&mut glossary_file, &glossary, SizeLimit::Infinite).unwrap();
            to_writer(&mut glossary_file, &glossary).unwrap();

            (glossary, keys)
        } else if let Some(glossary_filename) = options.value_of("use-glossary") {
            let key_filename = options.value_of("use-key");
            if key_filename.is_none() {
                println!("For now --use-glossary must go along with --use-key");
                return;
            }

            println!("Reading glossary");
            let mut glossary_file = File::open(glossary_filename).unwrap();
//             let glossary: Glossary = deserialize_from(&mut glossary_file, SizeLimit::Infinite).unwrap();
            let glossary: Glossary = from_reader(glossary_file).unwrap();

            println!("Reading key");
            let mut key_file = File::open(key_filename.unwrap()).unwrap();
//             let keys: KeyVec = deserialize_from(&mut key_file, SizeLimit::Infinite).unwrap();
            let keys: KeyVec = from_reader(&mut key_file).unwrap();

            (glossary, keys)
        } else {
            println!("Either build-glossary or use-glossary should be provided");
            return;
        }
    };

    if let Some(key_filename) = options.value_of("write-key") {
        println!("Writing key");

        let mut key_file = File::create(key_filename).unwrap();
//         serialize_into(&mut key_file, &keys, SizeLimit::Infinite).unwrap();
        to_writer(&mut key_file, &keys).unwrap();
    }

    if let Some(reconstruct_filename) = options.value_of("reconstruct") {
        println!("Reconstructing");

        sound::reconstruct(reconstruct_filename, &glossary, &keys, similarity);
    }
}
