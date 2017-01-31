
use std::collections::HashMap;

use std::collections::BTreeMap;
use std::collections::btree_map::Entry;
use std::collections::Bound::{Included, Excluded, Unbounded};

use std::rc::Rc;
use std::cell::RefCell;
use std::cmp::{Ord, Ordering, max};
use std::f32::consts::PI;

use std::fmt::{self, Debug, Formatter};

use serde::{Serialize, Deserialize, Serializer, Deserializer};
use serde::de::{Visitor, Error};

use sound; // {Spectrum, Cplx, Detector};
use sound::*;

// We need to implement serde for a foreign type
#[derive(Clone, Hash, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct BitVec(::bit_vec::BitVec);
use std::ops::{Deref, DerefMut};

impl BitVec {
    pub fn new() -> Self {
        BitVec(::bit_vec::BitVec::new())
    }

    pub fn from_bytes(bytes: &[u8]) -> Self {
        BitVec(::bit_vec::BitVec::from_bytes(bytes))
    }

    pub fn from_fn<F>(len: usize, f: F) -> Self
        where F: FnMut(usize) -> bool
    {
        BitVec(::bit_vec::BitVec::from_fn(len, f))
    }
}

impl Deref for BitVec {
    type Target = ::bit_vec::BitVec;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for BitVec {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Serialize for BitVec {
    fn serialize<S>(&self, serializer: &mut S) -> Result<(), S::Error>
        where S: Serializer
    {
        let bytes = self.to_bytes();
        serializer.serialize_bytes(&bytes)
    }
}

impl Deserialize for BitVec {
    fn deserialize<D>(deserializer: &mut D) -> Result<Self, D::Error>
        where D: Deserializer
    {
        struct BitVecVisitor;

        impl Visitor for BitVecVisitor {
            type Value = BitVec;

            #[inline]
            fn visit_bytes<E>(&mut self, v: &[u8]) -> Result<Self::Value, E>
                where E: Error
            {
                Ok(BitVec::from_bytes(v))
            }
        }

        deserializer.deserialize_bytes(BitVecVisitor)
    }
}


/// Special wrapper over `BitVec` that optimizes the case when
/// bit vector contains relatively small amount of set bits.
#[derive(Clone, Hash, Eq, PartialEq, Serialize, Deserialize)]
pub struct SparseBitVec {
    leading_zeros: usize,
    trailing_zeros: usize,
    bits_set: usize,
    bits: BitVec,
}

impl Debug for SparseBitVec {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{} {} {} [",
            self.leading_zeros,
            self.trailing_zeros,
            self.bits_set
        )?;

        for bit in self.bits.iter() {
            write!(f, "{}", if bit { "!" } else { "." } )?;
        }

        write!(f, "]")
    }
}

impl SparseBitVec {
    pub fn new() -> SparseBitVec {
        SparseBitVec {
            bits_set: 0,
            leading_zeros: 0,
            trailing_zeros: 0,
            bits: BitVec::new(),
        }
    }

    pub fn bits(&self) -> &BitVec {
        &self.bits
    }

    pub fn from_bitvec(bits: BitVec) -> SparseBitVec {
        let mut result = SparseBitVec {
            bits_set: 0,
            leading_zeros: 0,
            trailing_zeros: 0,
            bits: bits,
        };

        for bit in result.bits.iter() {
            if bit {
                result.bits_set += 1;
                result.trailing_zeros = 0;
            } else {
                if result.bits_set == 0 {
                    result.leading_zeros += 1;
                }
                result.trailing_zeros += 1;
            }
        }

        result
    }

    pub fn from_bytes(bytes: &[u8]) -> SparseBitVec {
        Self::from_bitvec(BitVec::from_bytes(bytes))
    }

    pub fn into_bitvec(self) -> BitVec {
        self.bits
    }

    pub fn fuzzy_eq(&self, other: &Self) -> usize {
        // Iterating through the vectors counting similarities
        // TODO Accept minimum similarity to speed up the merge

        let matched_bits = self
            .bits.blocks()
            .zip(other.bits.blocks())
            .fold(0, |count, (b1, b2)| count + (b1 & b2).count_ones() as usize);

        matched_bits
    }
}

impl PartialOrd for SparseBitVec {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Special implementation that is aware of vector's internal structure
impl Ord for SparseBitVec {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.leading_zeros < other.leading_zeros {
            Ordering::Greater
        } else if self.leading_zeros > other.leading_zeros {
            Ordering::Less
        } else {
            self.bits.cmp(&other.bits)
        }
    }
}

/// Sparse bit vector acting as a key of a fragment.
/// Type is used to differ fragment keys from other vectors.
#[derive(Clone, Hash, Eq, PartialEq, Ord, PartialOrd, Debug, Serialize, Deserialize)]
pub struct FragmentKey(SparseBitVec);

impl FragmentKey {
    pub fn new() -> FragmentKey {
        FragmentKey(SparseBitVec::new())
    }

    pub fn from_bitvec(bits: BitVec) -> FragmentKey {
        FragmentKey(SparseBitVec::from_bitvec(bits))
    }

    pub fn bits(&self) -> &BitVec {
        &self.0.bits()
    }
}

/// `Fragment` represents several time slices
/// of the spectrum within predefined frequency range.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Fragment {
    /// Slice of DFT within Dictionary's range.
    spectra: Vec<Spectrum>,

    /// Weight of fragment as prototype.
    /// Used during merge process.
    merge_weight: usize,

    // TODO use_count: Cell<usize>,
}

impl Fragment {
    pub fn new() -> Fragment {
        Fragment {
            spectra: Vec::new(),
            merge_weight: 1,
        }
    }

    pub fn from_spectra(spectra: Vec<Spectrum>) -> Fragment {
        Fragment {
            spectra: spectra,
            merge_weight: 1,
        }
    }

    pub fn spectrum(&self, slice_index: usize) -> &[Cplx] {
        &self.spectra[slice_index][..]
    }

    pub fn spectrum_mut(&mut self, slice_index: usize) -> &mut Spectrum {
        &mut self.spectra[slice_index]
    }

    pub fn spectra(&self) -> &Vec<Spectrum> {
        &self.spectra
    }
}

/// Internal container type used by `Dictionary`
type FragmentMap = BTreeMap<FragmentKey, Box<Fragment>>;

/// `Dictionary` holds fragments associated with bit vectors.
#[derive(Serialize, Deserialize)]
pub struct Dictionary {
    map: FragmentMap,
    lower_frequency: f32,
    upper_frequency: f32,
}

/// `Glossary` is a collection of `Dictionary`'s that spans
/// `the whole spectrum or it's interesting part.
#[derive(Serialize, Deserialize)]
pub struct Glossary {
    detectors: Vec<Detector>,
    dictionaries: Vec<Dictionary>,
}

impl Glossary {
    pub fn new(detectors: Vec<Detector>, dictionaries: Vec<Dictionary>) -> Glossary {
        Glossary {
            detectors: detectors,
            dictionaries: dictionaries,
        }
    }

    pub fn iter<'a>(&'a self) -> impl Iterator<Item=&'a Dictionary> {
        self.dictionaries.iter()
    }
}

impl Dictionary {
    pub fn len(&self) -> usize {
        self.map.len()
    }

    pub fn new(lower_frequency: f32, upper_frequency: f32) -> Dictionary {
        Dictionary {
            map: FragmentMap::new(),
            lower_frequency: lower_frequency,
            upper_frequency: upper_frequency,
        }
    }

    pub fn get_bounds(&self) -> (f32, f32) {
        (self.lower_frequency, self.upper_frequency)
    }

    pub fn lower_bound(key: &FragmentKey) -> FragmentKey {
        let mask = &key.0;
        let len = mask.bits.len();

        // Set all bits except trailing zeros to zero
        FragmentKey(SparseBitVec::from_bitvec(BitVec::from_fn(len, |x| {
            x > mask.leading_zeros + mask.bits_set
        })))
    }

    fn fragment_key(fragment_range: (f32, f32), detectors: &[Detector], fragment: &Fragment) -> FragmentKey {
        let mut result = BitVec::new();
        let base_index = (fragment_range.0 / BASE_FREQUENCY).round() as usize;

        // FIXME Should filter by detectors within the dictionary's range
        for spectrum in fragment.spectra.iter() {
            // Iterating through all detectors filtering out activity
            for detector in detectors {
                use sound::*;

                if detector.freq < fragment_range.0 || detector.freq > fragment_range.1 {
                    // Detector does not match frequency region of this dictionary
                    continue;
                }

                let detector_range = (fragment_range.0.max(detector.freq - detector.band), fragment_range.1.min(detector.freq + detector.band));

                // Each detector operates only in the fixed part of the spectrum
                // Selecting potentially interesting spectrum slice to check
                let low  = (detector_range.0.abs() / BASE_FREQUENCY).round() as usize - base_index;
                let high = (detector_range.1.abs() / BASE_FREQUENCY).round() as usize - base_index;

//                 println!("detector freq {} ± {}, amp {} ± {} dB, phase {}, spectrum len {}, {:?}",
//                     detector.freq,
//                     detector.band,
//                     detector.amp,
//                     AMPLITUDE_DEVIATION_DB,
//                     detector.phase,
//                     spectrum.len(),
//                     (low, high)
//                 );

                if low > spectrum.len() - 1 || high > spectrum.len() - 1 {
                    println!("invalid detector freq {}, band {}", detector.freq, detector.band);
                    break;
                }

                let range = &spectrum[low .. high + 1];

                // Selecting the entry with the largest amplitude
                // TODO vary sensitivity depending on the frequency deviation
                let (index, amplitude, phase) = range
                    .iter()
                    .enumerate()
                    .map(|(i, c)| (i, (c.norm() * 2.0) / NUM_POINTS as f32, c.arg()))
                    .max_by(|&(_, x, _), &(_, y, _)| float_cmp(x, y, 0.00001))
                    .unwrap_or((0, 0., 0.));

//                 let index = (detector.freq / BASE_FREQUENCY).round() as usize - base_index;
//                 let amplitude = (spectrum[index].norm() * 2.0) / NUM_POINTS as f32;
//                 let phase = spectrum[index].arg();

                let compress = |a| {
                    if a > -50. && a < -15. {
                        a + 10.
                    } else {
                        a
                    }
                };

                // Treating detector as active if max amplitude lays within detector's selectivity range
                let amp_match   = (compress(to_decibel(amplitude)).abs() - detector.amp.abs()).abs() < AMPLITUDE_DEVIATION_DB;

                // + PI is required to compare positive and negative values
                let phase_match = ((phase + PI) - (detector.phase + PI)).abs() < PHASE_DEVIATION_DB;

                let is_active   = amp_match && phase_match;
                result.push(is_active);

                if is_active {
//                     println!("[{}] detector freq {} ± {}, amp {} ± {} dB, phase {}",
//                         index,
//                         detector.freq,
//                         detector.band,
//                         detector.amp,
//                         AMPLITUDE_DEVIATION_DB,
//                         detector.phase
//                     );
//
//                     println!("signal amplitude {} → {} dB, phase {}{}\n",
//                         amplitude,
//                         to_decibel(amplitude),
//                         phase,
//                         if is_active { ", **MATCH**" } else { "" }
//                     );
                }
            }

            //sound::filter_detectors_inplace(&spectre, detectors, &mut result);
        }

        FragmentKey::from_bitvec(result)
    }

    /// Insert a fragment into the dictionary. If dictionary already
    /// contains fragment that is similar enough to the provided one
    /// then fragments are merged together. If no suitable match was
    /// found, then new item is inserted as is.
    pub fn insert_fragment(&mut self, fragment: Fragment, detectors: &[Detector], similarity: usize) -> Option<FragmentKey> {
        let freq_range = (self.lower_frequency, self.upper_frequency);
        let mut pending_key   = Self::fragment_key(freq_range, detectors, &fragment);
        let mut pending_value = Box::new(fragment);

        // Empty key means that the value may not be selected later
        if pending_key.0.bits_set == 0 {
            return None;
        }

        // Lower bound is the least meaningful element of the dictionary
        // which, if represented by a number, is less than the key's number
        let mut lower_bound = Self::lower_bound(&pending_key);

        enum KeyState {
            NotChanged,
            Changed(FragmentKey) // yields old key
        }

        loop {
            let key_state = {
                // Amount of bits to match the requested similarity percent
                let bit_threshold = (pending_key.0.bits_set as f32 / 100. * similarity as f32).round() as usize;

                // Finding best entry to merge-in the new value
                // TODO Optimize the case when at least threshold bits_set
                if let Some((key, value, _)) = self.map
                    .range_mut((Unbounded, Included(&pending_key)))
//                     .range_mut((Excluded(&lower_bound), Included(&pending_key)))
                    .map(|(k, v)| (k, v, k.0.fuzzy_eq(&pending_key.0)))
                    .filter(|&(_, _, m)| m >= bit_threshold)
                    .max_by(|x, y| x.2.cmp(&y.2)) // max by matched_bits
                {
                    // Best match is suitable for merge. Merging values
                    // and checking that the key wasn't changed during merge.
                    pending_key = Self::merge(freq_range, detectors, value, &pending_value);

                    // If key wasn't changed after merge then all is consistent
                    if *key == pending_key {
                        return Some(key.clone());
                    }

                    // Looks like key was changed after merge. We need to re-insert
                    // the merged value into the dictionary at the proper place
                    KeyState::Changed(key.clone())
                } else {
                    // Not enough matched bits to merge or dictionary is still empty
                    KeyState::NotChanged
                }
            };

            // If merge resulted in a changed key, then we need to reinsert
            // the value probably merging it again. Preparing for the next iteration.
            if let KeyState::Changed(old_key) = key_state {
                // Pending value will be re-inserted with a new pending key
                pending_value = self.map.remove(&old_key).unwrap();

                // Recalculating lower bound if new key has more bits to the right
                if pending_key.0.trailing_zeros < lower_bound.0.bits_set {
                    lower_bound = Self::lower_bound(&pending_key);
                }
            } else {
                // No suitable match was found, inserting fragment as the new prototype
                let key = pending_key.clone();
                self.map.insert(pending_key, pending_value);
                return Some(key);
            }
        }
    }

    pub fn find(&self, key: &FragmentKey, similarity: usize) -> Option<&Fragment> {
        if key.0.bits_set == 0 {
            return None;
        }

        // Lower bound is the least meaningful element of the dictionary
        // which, if represented by a number, is less than the key's number
        let lower_bound = Self::lower_bound(&key);

        // Amount of bits to match the requested similarity percent
        let bit_threshold = (key.0.bits_set as f32 / 100. * similarity as f32).round() as usize;

        if let Some((_, value, _)) = self.map
            .range((Unbounded, Included(key)))
//             .range((Excluded(&lower_bound), Included(key)))
            .map(|(k, v)| (k, v, k.0.fuzzy_eq(&key.0)))
            .filter(|&(_, _, m)| m >= bit_threshold)
            .max_by(|x, y| x.2.cmp(&y.2)) // max by matched_bits
        {
            Some(value)
        } else {
            None
        }
    }

    fn merge(freq_range: (f32, f32), detectors: &[Detector], prototype: &mut Fragment, value: &Fragment) -> FragmentKey {
        assert_eq!(prototype.spectra.len(), value.spectra.len());

        let dest_weight = prototype.merge_weight as f32 / (prototype.merge_weight + value.merge_weight) as f32;
        let src_weight  = value.merge_weight as f32 / (prototype.merge_weight + value.merge_weight) as f32;

        // For each spectrum slice merging complex values by weight
        for (dest, src) in prototype.spectra.iter_mut().zip(value.spectra.iter()) {
            assert_eq!(dest.len(), src.len());

            for (c1, c2) in dest.iter_mut().zip(src.iter()) {
                let new_norm  = c1.norm() * dest_weight + c2.norm() * src_weight;
                let new_tetha = c1.arg()  * dest_weight + c2.arg()  * src_weight;

                *c1 = Cplx::from_polar(&new_norm, &new_tetha);
            }
        }

        // Prototype now accumulates both values, so it's weight increases
        prototype.merge_weight += value.merge_weight;

        // TODO Calculate during merge
        Self::fragment_key(freq_range, &detectors, &prototype)
    }

    pub fn iter<'b>(&'b self) -> impl Iterator<Item=(&'b FragmentKey, &'b Box<Fragment>)> {
        self.map.iter()
    }
}

#[cfg(test)]
mod sparse_bitvec {
    use super::*;

    #[test] fn new() {
        let vec = SparseBitVec::new();
        assert_eq!(vec.bits_set, 0);
        assert_eq!(vec.leading_zeros, 0);
        assert_eq!(vec.trailing_zeros, 0);
    }

    #[test] fn from_bitvec() {
        let vec = SparseBitVec::from_bytes(&[
            0b_0000_0001,
            0b_1010_0000,
            0b_0001_0011,
            0b_0000_0000
        ]);

        assert_eq!(vec.bits_set, 6);
        assert_eq!(vec.leading_zeros, 7);
        assert_eq!(vec.trailing_zeros, 8);
    }

    #[test] fn zeros() {
        let data = &[
            // mask          leading  trailing
            ([0b_1111_1111], 0,       0),
            ([0b_1000_0001], 0,       0),
            ([0b_1000_0000], 0,       7),
            ([0b_0000_0001], 7,       0),
            ([0b_0100_0010], 1,       1),
            ([0b_0010_0100], 2,       2),
            ([0b_0001_1000], 3,       3),
            ([0b_0000_0000], 8,       8),
        ];

        for &(mask, leading, trailing) in data {
            let vec = SparseBitVec::from_bytes(&mask);

            assert_eq!(vec.leading_zeros, leading);
            assert_eq!(vec.trailing_zeros, trailing);
        }
    }

    #[test] fn set_bits() {
        let data = &[
            // mask         bits set
            ([0b_0000_0000], 0),
            ([0b_1111_1111], 8),
            ([0b_1010_1010], 4),
            ([0b_0101_0101], 4),

            ([0b_1000_0000], 1),
            ([0b_0000_0001], 1),

            ([0b_1000_0001], 2),
            ([0b_0001_1000], 2),

            ([0b_1000_0000], 1),
            ([0b_1000_0000], 1),

        ];

        for &(mask, bits) in data {
            let vec = SparseBitVec::from_bytes(&mask);

            assert_eq!(vec.bits_set, bits);
        }
    }

    #[test] fn cmp() {
        let data = &[
            // first vector  second vector   order
            ([0b_0000_0000], [0b_0000_0000], Ordering::Equal),
            ([0b_1111_1111], [0b_1111_1111], Ordering::Equal),
            ([0b_0111_1111], [0b_0111_1111], Ordering::Equal),
            ([0b_1111_1110], [0b_1111_1110], Ordering::Equal),
            ([0b_0001_1000], [0b_0001_1000], Ordering::Equal),

            ([0b_0000_0000], [0b_0000_0001], Ordering::Less),
            ([0b_1111_1110], [0b_1111_1111], Ordering::Less),
            ([0b_1111_0000], [0b_1111_0001], Ordering::Less),
            ([0b_1111_1110], [0b_1111_1111], Ordering::Less),
            ([0b_0011_1100], [0b_0011_1110], Ordering::Less),
            ([0b_0101_0101], [0b_1010_1010], Ordering::Less),

            ([0b_0000_0001], [0b_0000_0000], Ordering::Greater),
            ([0b_1111_1111], [0b_1111_1110], Ordering::Greater),
            ([0b_1111_0001], [0b_1111_0000], Ordering::Greater),
            ([0b_1111_1111], [0b_1111_1110], Ordering::Greater),
            ([0b_0011_1110], [0b_0011_1100], Ordering::Greater),
            ([0b_1010_1010], [0b_0101_0101], Ordering::Greater),
        ];

        for &(b1, b2, order) in data {
            let v1 = SparseBitVec::from_bytes(&b1);
            let v2 = SparseBitVec::from_bytes(&b2);

            assert_eq!(v1.cmp(&v2), order, "{:?} vs {:?}", v1, v2);
            assert_eq!(v1.partial_cmp(&v2), Some(order), "{:?} vs {:?}", v1, v2);
        }
    }

    #[test] fn cmp_empty() {
        assert_eq!(SparseBitVec::new(), SparseBitVec::new());
    }
}
