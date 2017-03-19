
#![feature(btree_range, collections_bound)]
#![feature(conservative_impl_trait)]
#![allow(dead_code)]

extern crate num_traits;
extern crate num_complex;
extern crate bit_vec;
extern crate dft;
extern crate clap;
extern crate hound;
extern crate sample;
extern crate rayon;
extern crate itertools;

#[macro_use] extern crate serde_derive;
extern crate serde;

extern crate bincode;
use bincode::{serialize_into, deserialize_from};

use clap::{Arg, App};

mod sound;
mod memory;
mod iter;

use sound::*;
use memory::*;

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

    let similarity = options.value_of("similarity").unwrap_or("50").parse().unwrap();

    let (glossary, keys) = {
        if let Some(glossary_filename) = options.value_of("build-glossary") {
            let input_filename = options.value_of("input").unwrap();
            let glossary = sound::build_glossary(input_filename, similarity);

            println!("Writing glossary file");
            let mut glossary_file = File::create(glossary_filename).unwrap();
            serialize_into(&mut glossary_file, &glossary, bincode::Infinite).unwrap();
            // to_writer(&mut glossary_file, &glossary).unwrap();


            let keys = if let Some(key_filename) = options.value_of("write-key") {
                let keys = sound::analyze_file(input_filename);

                println!("Writing key");
                let mut key_file = File::create(key_filename).unwrap();
                serialize_into(&mut key_file, &keys, bincode::Infinite).unwrap();
                // to_writer(&mut key_file, &keys).unwrap();

                Some(keys)
            } else {
                None
            };

            (glossary, keys)
        } else if let Some(glossary_filename) = options.value_of("use-glossary") {
            let key_filename = options.value_of("use-key");
            if key_filename.is_none() {
                println!("For now --use-glossary must go along with --use-key");
                return;
            }

            println!("Reading glossary");
            let mut glossary_file = File::open(glossary_filename).unwrap();
            let glossary: Glossary = deserialize_from(&mut glossary_file, bincode::Infinite).unwrap();
            // let glossary: Glossary = from_reader(glossary_file).unwrap();

            println!("Reading key");
            let mut key_file = File::open(key_filename.unwrap()).unwrap();
            let keys: KeyVec = deserialize_from(&mut key_file, bincode::Infinite).unwrap();
            // let keys: KeyVec = from_reader(&mut key_file).unwrap();

            (glossary, Some(keys))
        } else {
            println!("Either build-glossary or use-glossary should be provided");
            return;
        }
    };

    if let Some(reconstruct_filename) = options.value_of("reconstruct") {
        println!("Reconstructing");

        sound::reconstruct(reconstruct_filename, &glossary, &keys.unwrap(), 5);
    }
}
