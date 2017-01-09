use std::collections::HashMap;
use std::collections::hash_map::Entry;

use std::rc::Rc;
use std::cell::RefCell;

use bit_vec::BitVec;

use sound::{Spectrum, Cplx, Detector};

#[derive(Hash, Eq, PartialEq, Debug)]
struct FragmentKey(BitVec);

const FRAGMENT_SLICES: usize = 8;

/// `Fragment` represents several time slices
/// of the spectrum within specified range.
struct Fragment {
    spectra: Vec<Spectrum>,
    merge_factor: i32,
}

/// Internal container type used by `Dictionary`
type RcFragment  = Rc<Fragment>;
type FragmentMap = HashMap<FragmentKey, RcFragment>;

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
    fn insert(&mut self, fragment: Fragment) {
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
    }

    fn find(&self, key: &FragmentKey) -> Option<RcFragment> {
        match self.map.get(key) {
            Some(fragment) => Some(fragment.clone()),
            None => None
        }
    }

    fn merge(prototype: &mut Fragment, new_prototype: &Fragment) -> FragmentKey {
        prototype.merge_factor += new_prototype.merge_factor;

        FragmentKey::new()
    }
}
