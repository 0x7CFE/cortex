use std::collections::BTreeMap;
use std::collections::btree_map::Entry;

use std::rc::Rc;
use std::cell::RefCell;
use std::cmp::{Ord, Ordering};

use bit_vec::BitVec;

use sound::{Spectrum, Cplx, Detector};


#[derive(Hash, Eq, PartialEq, Debug)]
struct SparseBitVec {
    leading_zeros: usize,
    trailing_zeros: usize,
    bits_set: usize,
    data: BitVec,
}

impl SparseBitVec {
    fn new() -> SparseBitVec {
        SparseBitVec {
            bits_set: 0,
            leading_zeros: 0,
            trailing_zeros: 0,
            data: BitVec::new(),
        }
    }

    fn bits(&self) -> &BitVec {
        &self.data
    }

    fn from_bitvec(bits: BitVec) -> SparseBitVec {
        let mut idea = SparseBitVec {
            bits_set: 0,
            leading_zeros: 0,
            trailing_zeros: 0,
            data: bits,
        };

        for bit in idea.data.iter() {
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
        self.data
    }
}

impl PartialOrd for SparseBitVec {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SparseBitVec {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.leading_zeros < other.leading_zeros {
            Ordering::Less
        } else if self.leading_zeros > other.leading_zeros {
            Ordering::Greater
        } else {
            self.data.cmp(&other.data)
        }
    }
}

#[derive(Hash, Eq, PartialEq, Ord, PartialOrd, Debug)]
struct FragmentKey(SparseBitVec);

/// `Fragment` represents several time slices
/// of the spectrum within specified range.
struct Fragment {
    spectra: Vec<Spectrum>,
    merge_factor: i32,
}

const FRAGMENT_SLICES: usize = 8;

/// Internal container type used by `Dictionary`
type RcFragment  = Rc<Fragment>;
type FragmentMap = BTreeMap<FragmentKey, Box<Fragment>>;

/// `Dictionary` holds fragments associated with bit vectors.
struct Dictionary<'a> {
    map: FragmentMap,
    detectors: &'a [Detector]
}

impl FragmentKey {
    fn new() -> FragmentKey {
        FragmentKey(SparseBitVec::new())
    }

    fn from_bitvec(bits: BitVec) -> FragmentKey {
        FragmentKey(SparseBitVec::from_bitvec(bits))
    }
}

impl Fragment {
    fn new() -> Fragment {
        Fragment {
            spectra: Vec::new(),
            merge_factor: 1,
        }
    }

    fn spectrum(&self, slice_index: usize) -> &[Cplx] {
        &self.spectra[slice_index][..]
    }

    fn spectrum_mut(&mut self, slice_index: usize) -> &mut Spectrum {
        &mut self.spectra[slice_index]
    }

    fn key(&self, detectors: &[Detector]) -> FragmentKey {
        FragmentKey::new()
    }
}

impl<'a> Dictionary<'a> {
    fn new(detectors: &[Detector]) -> Dictionary {
        Dictionary {
            map: FragmentMap::new(),
            detectors: detectors
        }
    }

    /// Insert a fragment into the dictionary. If dictionary already
    /// contains fragment that is similar enough to the one provided
    /// then fragments are merged together and the hash of result is
    /// returned. If no match was found, then new item is added.
    /*fn insert(&mut self, fragment: Fragment) {
        let mut new_key = fragment.key(self.detectors);
        let mut new_prototype = Rc::new(fragment);

        loop {
            match self.map.entry(new_key) {
                Entry::Vacant(entry) => {
                    // Fragment is unique, inserting it as prototype
                    entry.insert(new_prototype);
                    break;
                },

                Entry::Occupied(mut entry) => {
                    let key_differs = {
                        new_key = Self::merge(Rc::get_mut(entry.get_mut()).unwrap(), &new_prototype);
                        *entry.key() != new_key
                    };

                    if key_differs {
                        // Need to re-insert with a new key
                        let (_, prototype) = entry.remove_entry();
                        new_prototype = prototype;
                    } else {
                        // Key is the same, all is OK
                        break;
                    }
                }
            }
        }
    }*/

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

    fn merge(prototype: &mut Fragment, new_prototype: &Fragment) -> FragmentKey {
        prototype.merge_factor += new_prototype.merge_factor;

        FragmentKey::new()
    }
}
