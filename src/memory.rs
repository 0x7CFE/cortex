
use std::collections::BTreeMap;
use std::collections::btree_map::Entry;
use std::collections::Bound::{Included, Excluded, Unbounded};

use std::rc::Rc;
use std::cell::RefCell;
use std::cmp::{Ord, Ordering, max};

use bit_vec::BitVec;

use sound::{Spectrum, Cplx, Detector};

/// Special wrapper over `BitVec` that optimizes the case when
/// bit vector contains relatively small amount of set bits.
#[derive(Clone, Hash, Eq, PartialEq, Debug)]
pub struct SparseBitVec {
    leading_zeros: usize,
    trailing_zeros: usize,
    bits_set: usize,
    bits: BitVec,
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
        let mut idea = SparseBitVec {
            bits_set: 0,
            leading_zeros: 0,
            trailing_zeros: 0,
            bits: bits,
        };

        for bit in &idea.bits {
            if bit {
                idea.bits_set += 1;
                idea.trailing_zeros = 0;
            } else {
                if idea.bits_set == 0 {
                    idea.leading_zeros += 1;
                }
                idea.trailing_zeros += 1;
            }
        }

        idea
    }

    pub fn from_bytes(bytes: &[u8]) -> SparseBitVec {
        Self::from_bitvec(BitVec::from_bytes(bytes))
    }

    pub fn into_bitvec(self) -> BitVec {
        self.bits
    }

    pub fn fuzzy_eq(&self, other: &Self) -> usize {
        // Iterating through the vectors counting similarities
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
#[derive(Clone, Hash, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub struct FragmentKey(SparseBitVec);

impl FragmentKey {
    pub fn new() -> FragmentKey {
        FragmentKey(SparseBitVec::new())
    }

    pub fn from_bitvec(bits: BitVec) -> FragmentKey {
        FragmentKey(SparseBitVec::from_bitvec(bits))
    }
}

/// Amount of spectrum slices per single fragment
pub const SPECTRA_PER_FRAGMENT: usize = 8;

/// `Fragment` represents several time slices
/// of the spectrum within specified range.
pub struct Fragment {
    /// Complete slice of DFT.
    spectra: Vec<Spectrum>,

    /// Weight of fragment as prototype.
    /// Used during merge process.
    merge_weight: i32,
}

/// Internal container type used by `Dictionary`
type FragmentMap = BTreeMap<FragmentKey, Box<Fragment>>;

/// `Dictionary` holds fragments associated with bit vectors.
pub struct Dictionary<'a> {
    map: FragmentMap,
    detectors: &'a [Detector]
}

impl Fragment {
    pub fn new() -> Fragment {
        Fragment {
            spectra: Vec::with_capacity(SPECTRA_PER_FRAGMENT),
            merge_weight: 1,
        }
    }

    pub fn spectrum(&self, slice_index: usize) -> &[Cplx] {
        &self.spectra[slice_index][..]
    }

    pub fn spectrum_mut(&mut self, slice_index: usize) -> &mut Spectrum {
        &mut self.spectra[slice_index]
    }

    pub fn key(&self, detectors: &[Detector]) -> FragmentKey {
        let mut result = BitVec::new();

        for spectre in self.spectra.iter() {
            sound::filter_detectors_inplace(&spectre, detectors, &mut result);
        }

        FragmentKey::from_bitvec(result)
    }
}

impl<'a> Dictionary<'a> {
    pub fn new(detectors: &[Detector]) -> Dictionary {
        Dictionary {
            map: FragmentMap::new(),
            detectors: detectors
        }
    }

    pub fn lower_bound(key: &FragmentKey) -> FragmentKey {
        let mask = &key.0;
        let len = mask.bits.len();

        // Set all bits except trailing zeros to zero
        FragmentKey(SparseBitVec::from_bitvec(BitVec::from_fn(len, |x| {
            x > mask.leading_zeros + mask.bits_set
        })))
    }

    /// Insert a fragment into the dictionary. If dictionary already
    /// contains fragment that is similar enough to the provided one
    /// then fragments are merged together. If no suitable match was
    /// found, then new item is inserted as is.
    pub fn insert(&mut self, fragment: Fragment, similarity: usize) {

        let mut pending_key = fragment.key(self.detectors);
        let mut pending_value = Box::new(fragment);

        // FIXME Should we treat is as NOOP if empty key is pending?
        assert!(pending_key.0.bits_set > 0);

        // Lower bound is the least meaningful element of the dictionary
        // which, if represented by a number, is less than the key's number
        let mut lower_bound = Self::lower_bound(&pending_key);

        loop {
            let key_changed = {
                // Finding best entry to merge-in the new value
                if let Some((ref key, ref mut value, matched_bits)) = self.map
                    .range_mut(Excluded(&lower_bound), Included(&pending_key))
                    .map(|(k, v)| (k, v, k.0.fuzzy_eq(&pending_key.0)))
                    .filter(|&(_, _, m)| m >= similarity)
                    .max_by(|x, y| x.2.cmp(&y.2)) // max by matched_bits
                {
                    // Best match is suitable for merge. Merging values
                    // and checking that the key wasn't changed during merge.
                    pending_key = Self::merge(value, &pending_value);

                    // If key wasn't changed after merge then all is consistent
                    if **key == pending_key {
                        return;
                    }

                    // Looks like key was changed after merge. We need to re-insert
                    // the merged value into the dictionary at the proper place
                    true
                } else {
                    // Not enough matched bits to merge or dictionary is still empty
                    false
                }
            };

            // If merge resulted in a changed key, then we need to reinsert
            // the value probably merging it again. Preparing for the next iteration.
            if key_changed {
                // New key now stores key after merge
                pending_value = self.map.remove(&pending_key).unwrap();

                // Recalculating lower bound if new key has more bits to the right
                if pending_key.0.trailing_zeros < lower_bound.0.bits_set {
                    lower_bound = Self::lower_bound(&pending_key);
                }
            } else {
                // No suitable match was found, inserting fragment as the new prototype
                self.map.insert(pending_key, pending_value);
                return;
            }
        }
    }

    pub fn find(&self, key: &FragmentKey, similarity: usize) -> Option<&Fragment> {
        // Lower bound is the least meaningful element of the dictionary
        // which, if represented by a number, is less than the key's number
        let lower_bound = Self::lower_bound(&key);

        if let Some((_, ref value, _)) = self.map
            .range(Excluded(&lower_bound), Included(&key))
            .map(|(k, v)| (k, v, k.0.fuzzy_eq(&key.0)))
            .filter(|&(_, _, m)| m >= similarity)
            .max_by(|x, y| x.2.cmp(&y.2)) // max by matched_bits
        {
            Some(&value)
        } else {
            None
        }
    }

    fn merge(prototype: &mut Fragment, pending_value: &Fragment) -> FragmentKey {
        prototype.merge_weight += pending_value.merge_weight;

        FragmentKey::new() // TODO
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
        let plan = vec![
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

        for (mask, leading, trailing) in plan {
            let vec = SparseBitVec::from_bytes(&mask);

            assert_eq!(vec.leading_zeros, leading);
            assert_eq!(vec.trailing_zeros, trailing);
        }
    }

    #[test] fn set_bits() {
        let plan = vec![
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

        for (mask, bits) in plan {
            let vec = SparseBitVec::from_bytes(&mask);

            assert_eq!(vec.bits_set, bits);
        }
    }

    #[test] fn cmp() {
        let plan = vec![
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

        for (b1, b2, order) in plan {
            let v1 = SparseBitVec::from_bytes(&b1);
            let v2 = SparseBitVec::from_bytes(&b2);

            assert_eq!(v1.cmp(&v2), order, "{:?} vs {:?}", v1, v2);
            assert_eq!(v1.partial_cmp(&v2), Some(order), "{:?} vs {:?}", v1, v2);
        }
    }
}
}
