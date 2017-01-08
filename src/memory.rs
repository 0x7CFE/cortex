use sound::{Spectrum, Cplx, Detector};
use bit_vec::BitVec;
use std::collections::HashMap;
use std::collections::hash_map::Entry;

#[derive(Hash, Eq, PartialEq, Debug)]
struct FragmentKey(BitVec);

/// `Fragment` represents several time slices
/// of the spectrum within specified range.
struct Fragment {
    merge_factor: i32,
    spectra: Vec<Spectrum>,
}

/// Internal container type used by `Dictionary`
type FragmentMap = HashMap<FragmentKey, Fragment>;

/// `Dictionary` holds fragments associated with bit vectors.
struct Dictionary<'a> {
    map: FragmentMap,
    detectors: &'a [Detector]
}

impl FragmentKey {
    fn new() -> FragmentKey {
        FragmentKey(BitVec::new())
    }
}

impl Fragment {
    fn new() -> Fragment {
        Fragment {
            merge_factor: 1,
            spectra: Vec::new()
        }
    }

    fn slice(&self, index: usize) -> &[Cplx] {
        &self.spectra[index][..]
    }

    fn slice_mut(&mut self, index: usize) -> &mut Spectrum {
        &mut self.spectra[index]
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
    fn insert(&mut self, fragment: Fragment) {
        let mut new_key = fragment.key(self.detectors);
        let mut new_prototype = fragment;

        loop {
            match self.map.entry(new_key) {
                Entry::Vacant(entry) => {
                    // Fragment is unique, inserting it as prototype
                    entry.insert(new_prototype);
                    break;
                },

                Entry::Occupied(mut entry) => {
                    let key_differs = {
                        new_key = Self::merge(entry.get_mut(), &new_prototype);
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
    }

    fn find(&self, key: &FragmentKey) -> Option<&Fragment> {
        self.map.get(key)
    }

    fn merge(prototype: &mut Fragment, new_prototype: &Fragment) -> FragmentKey {
        prototype.merge_factor += new_prototype.merge_factor;

        FragmentKey::new()
    }
}
