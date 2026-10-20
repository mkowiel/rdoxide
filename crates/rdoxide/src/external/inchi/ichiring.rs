//! Ring detection by BFS (`INCHI_BASE/src/ichiring.{h,c}`): the queue and the
//! smallest ring through a bond.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use super::inpdef::{InpAtom, MAX_ATOMS};

/// `QUEUE` of `qInt` (`AT_RANK`) values: a ring buffer of fixed size.
#[derive(Debug, Clone)]
pub struct Queue {
    val: Vec<u16>,
    n_tot_length: i32,
    /// element to remove if `n_length > 0`
    n_first: i32,
    /// `n_first + n_length` is the next free position
    n_length: i32,
}

impl Queue {
    /// `QueueCreate`: `None` for a length below 1.
    pub fn create(n_tot_length: i32) -> Option<Queue> {
        if n_tot_length < 1 {
            return None;
        }
        Some(Queue {
            val: vec![0; n_tot_length as usize],
            n_tot_length,
            n_first: 0,
            n_length: 0,
        })
    }

    /// `QueueAdd`: the new length, or -1 if the queue is full.
    pub fn add(&mut self, val: u16) -> i32 {
        if self.n_length < self.n_tot_length {
            let pos = (self.n_first + self.n_length) % self.n_tot_length;
            self.val[pos as usize] = val;
            self.n_length += 1;
            return self.n_length;
        }
        -1
    }

    /// `QueueGet`: reads and removes the first value; returns the remaining
    /// length, or -1 if the queue is empty.
    pub fn get(&mut self) -> Result<(u16, i32), i32> {
        if self.n_length > 0 {
            let v = self.val[self.n_first as usize];
            // the retrieved value is not overwritten
            self.n_first = if self.n_first == self.n_tot_length - 1 {
                0
            } else {
                self.n_first + 1
            };
            self.n_length -= 1;
            return Ok((v, self.n_length));
        }
        Err(-1)
    }

    /// `QueueGetAny`: the value at the internal position `ord`.
    pub fn get_any(&self, ord: i32) -> Option<u16> {
        if 0 <= ord && ord < self.n_tot_length {
            Some(self.val[ord as usize])
        } else {
            None
        }
    }

    /// `QueueReinit`: empties the queue; returns its total length.
    pub fn reinit(&mut self) -> i32 {
        self.n_first = 0;
        self.n_length = 0;
        self.n_tot_length
    }

    /// `QueueLength`
    pub fn length(&self) -> i32 {
        self.n_length
    }

    /// `QueueWrittenLength`: the number of used internal elements.
    pub fn written_length(&self) -> i32 {
        let len = self.n_first + self.n_length;
        len.min(self.n_tot_length)
    }
}

/// `GetMinRingSize`: BFS from the queued atoms; the smallest ring found, 0 if
/// none is smaller than `n_max_ring_size`, or -1 on an error.
fn get_min_ring_size(
    atom: &[InpAtom],
    q: &mut Queue,
    n_atom_level: &mut [u16],
    c_source: &mut [i8],
    n_max_ring_size: u16,
) -> i32 {
    let max_atoms = MAX_ATOMS as u16;
    let mut n_min_ring_size: u16 = max_atoms + 1;
    loop {
        let q_len = q.length();
        if q_len == 0 {
            break;
        }
        // traverse the next level (the next outer ring)
        for _ in 0..q_len {
            let Ok((at_no, _)) = q.get() else {
                return -1;
            };
            let iat_no = at_no as usize;
            let n_cur_level = n_atom_level[iat_no] + 1;
            if 2 * u32::from(n_cur_level) > u32::from(n_max_ring_size) + 4 {
                // 2 * level = ring size + 3 + k, k = 0 or 1
                if n_min_ring_size < max_atoms + 1 {
                    return if n_min_ring_size >= n_max_ring_size {
                        0
                    } else {
                        i32::from(n_min_ring_size)
                    };
                }
                // the smallest ring is larger than n_max_ring_size
                return 0;
            }
            for j in 0..atom[iat_no].valence as usize {
                let next = atom[iat_no].neighbor[j];
                let inext = next as usize;
                if n_atom_level[inext] == 0 {
                    // not traversed yet: add to the queue
                    if q.add(next) >= 0 {
                        n_atom_level[inext] = n_cur_level;
                        // keep the path number
                        c_source[inext] = c_source[iat_no];
                    } else {
                        return -1;
                    }
                } else if u32::from(n_atom_level[inext]) + 1 >= u32::from(n_cur_level)
                    && c_source[inext] != c_source[iat_no]
                {
                    // a ring closure
                    if c_source[inext] == -1 {
                        return -1;
                    }
                    let n_ring_size = n_atom_level[inext]
                        .wrapping_add(n_cur_level)
                        .wrapping_sub(2);
                    if n_ring_size < n_min_ring_size {
                        n_min_ring_size = n_ring_size;
                    }
                }
            }
        }
    }
    if n_min_ring_size < max_atoms + 1 {
        return if n_min_ring_size >= n_max_ring_size {
            0
        } else {
            i32::from(n_min_ring_size)
        };
    }
    0
}

/// `is_bond_in_Nmax_memb_ring`: the size of the smallest ring through the
/// bond `atom[at_no].neighbor[neigh_ord]` if it is smaller than
/// `n_max_ring_size`; 0 if `n_max_ring_size < 3` or no such ring; negative on
/// an error. `n_atom_level` and `c_source` are work arrays (all zero on
/// entry and on return).
pub fn is_bond_in_nmax_memb_ring(
    atom: &[InpAtom],
    at_no: usize,
    neigh_ord: usize,
    q: &mut Queue,
    n_atom_level: &mut [u16],
    c_source: &mut [i8],
    n_max_ring_size: u16,
) -> i32 {
    if n_max_ring_size < 3 {
        return 0;
    }
    q.reinit();
    // mark the starting atom
    n_atom_level[at_no] = 1;
    c_source[at_no] = -1;
    // add the neighbors
    for i in 0..atom[at_no].valence as usize {
        let n = atom[at_no].neighbor[i];
        n_atom_level[n as usize] = 2;
        c_source[n as usize] = 1 + i8::from(i == neigh_ord);
        q.add(n);
    }
    let n_min_ring_size = get_min_ring_size(atom, q, n_atom_level, c_source, n_max_ring_size);
    // cleanup
    for i in 0..q.written_length() {
        if let Some(n) = q.get_any(i) {
            n_atom_level[n as usize] = 0;
            c_source[n as usize] = 0;
        }
    }
    n_atom_level[at_no] = 0;
    c_source[at_no] = 0;
    n_min_ring_size
}

/// `is_atom_in_3memb_ring`
pub fn is_atom_in_3memb_ring(atom: &[InpAtom], at_no: usize) -> bool {
    if atom[at_no].n_num_at_in_ring_system < 3 {
        return false;
    }
    let val = atom[at_no].valence as usize;
    for i in 0..val {
        let neigh = atom[at_no].neighbor[i] as usize;
        if atom[at_no].n_ring_system != atom[neigh].n_ring_system {
            continue;
        }
        for j in 0..atom[neigh].valence as usize {
            let neigh_neigh = atom[neigh].neighbor[j];
            if neigh_neigh as usize == at_no {
                continue;
            }
            if atom[at_no].neighbor[..val].contains(&neigh_neigh) {
                return true;
            }
        }
    }
    false
}
