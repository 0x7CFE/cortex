
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
struct SparseBitVec {
    leading_zeros: usize,
    trailing_zeros: usize,
    bits_set: usize,
    bits: BitVec,
}

impl SparseBitVec {
    fn new() -> SparseBitVec {
        SparseBitVec {
            bits_set: 0,
            leading_zeros: 0,
            trailing_zeros: 0,
            bits: BitVec::new(),
        }
    }

    fn bits(&self) -> &BitVec {
        &self.bits
    }

    fn from_bitvec(bits: BitVec) -> SparseBitVec {
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

    fn into_bitvec(self) -> BitVec {
        self.bits
    }

    fn fuzzy_eq(&self, other: &Self) -> usize {
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
struct FragmentKey(SparseBitVec);

impl FragmentKey {
    fn new() -> FragmentKey {
        FragmentKey(SparseBitVec::new())
    }

    fn from_bitvec(bits: BitVec) -> FragmentKey {
        FragmentKey(SparseBitVec::from_bitvec(bits))
    }
}

/// Amount of spectrum slices per single fragment
const SPECTRA_PER_FRAGMENT: usize = 8;

/// `Fragment` represents several time slices
/// of the spectrum within specified range.
struct Fragment {
    /// Complete slice of DFT.
    spectra: Vec<Spectrum>,

    /// Weight of fragment as prototype.
    /// Used during merge process.
    merge_weight: i32,
}

/// Internal container type used by `Dictionary`
type FragmentMap = BTreeMap<FragmentKey, Box<Fragment>>;

/// `Dictionary` holds fragments associated with bit vectors.
struct Dictionary<'a> {
    map: FragmentMap,
    detectors: &'a [Detector]
}

impl Fragment {
    fn new() -> Fragment {
        Fragment {
            spectra: Vec::with_capacity(SPECTRA_PER_FRAGMENT),
            merge_weight: 1,
        }
    }

    fn spectrum(&self, slice_index: usize) -> &[Cplx] {
        &self.spectra[slice_index][..]
    }

    fn spectrum_mut(&mut self, slice_index: usize) -> &mut Spectrum {
        &mut self.spectra[slice_index]
    }

    fn key(&self, detectors: &[Detector]) -> FragmentKey {
        FragmentKey::new() // TODO
    }
}

impl<'a> Dictionary<'a> {
    fn new(detectors: &[Detector]) -> Dictionary {
        Dictionary {
            map: FragmentMap::new(),
            detectors: detectors
        }
    }

    fn lower_bound(key: &FragmentKey) -> FragmentKey {
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
    fn insert(&mut self, fragment: Fragment, similarity: usize) {

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
                if let Some((key, value, matched_bits)) = self.map
                        .range_mut(Excluded(&lower_bound), Included(&pending_key))
                        .map(|(k, v)| (k, v, k.0.fuzzy_eq(&pending_key.0)))
                        .filter(|&(_, _, m)| m >= similarity)
                        .max_by(|x, y| x.2.cmp(&y.2)) // max by matched_bits
                {
                    // Best match is suitable for merge. Merging values
                    // and checking that the key wasn't changed during merge.
                    pending_key = Self::merge(value, &pending_value);

                    // If key wasn't changed after merge then all is consistent
                    if *key == pending_key {
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

    fn find(&self, key: &FragmentKey, similarity: u32) -> Option<&Fragment> {
        //let mask = &key.0.bits();


        match self.map.get(key) {
            Some(fragment) => Some(&fragment),
            None => None
        }

        /*for &(FragmentKey(ref candidate), ref value) in self.map.iter() {
            assert_eq!(mask.len(), candidate.bits().len());

            // Iterating through the keys looking for similarities.
            // When there are more than `similarity` matched bits
            // keys are said to be matching. Otherwise it is a no-match.
            let mut total_match = 0;
            for (b1, b2) in mask.blocks().zip(candidate.bits().blocks()) {
                total_match += (b1 & b2).count_ones();

                if total_match >= similarity {
                    return Some(value);
                }
            }
        }

        None*/
    }

    fn merge(prototype: &mut Fragment, pending_value: &Fragment) -> FragmentKey {
        prototype.merge_weight += pending_value.merge_weight;

        FragmentKey::new() // TODO
    }
}
