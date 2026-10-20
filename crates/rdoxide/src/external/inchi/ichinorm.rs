//! Ring systems and component marks of `inp_ATOM` arrays
//! (`INCHI_BASE/src/ichinorm.c`: `MarkRingSystemsInp`, `UnMark*`;
//! `UnMarkRingSystemsInp` is in `runichi3.c`).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).
//!
//! `FIND_RINS_SYSTEMS_DISTANCES` is 0 in the library build, so the distance
//! from terminal ring systems is not computed. The underivatization and
//! ring-chain code of the file (engineering mode only) is not ported.

use super::ichierr::CT_OVERFLOW;
use super::inpdef::InpAtom;
use super::strutil::OrigAtomData;

/// `subtract_DT_from_num_H`: removes the isotopic H counts from `num_H`
pub fn subtract_dt_from_num_h(at: &mut [InpAtom]) -> i32 {
    for a in at.iter_mut() {
        for j in 0..a.num_iso_h.len() {
            a.num_h -= a.num_iso_h[j];
        }
    }
    0
}

/// `UnMarkDisconnectedComponents`
pub fn unmark_disconnected_components(orig: &mut OrigAtomData) -> i32 {
    for a in orig.at.iter_mut().take(orig.num_inp_atoms) {
        a.orig_compt_at_numb = 0;
        a.component = 0;
    }
    orig.n_cur_at_len = None;
    orig.n_old_comp_number = None;
    orig.num_components = 0;
    0
}

/// `UnMarkOtherIndicators`
pub fn unmark_other_indicators(at: &mut [InpAtom]) -> i32 {
    for a in at.iter_mut() {
        a.at_type = 0;
        a.c_flags = 0;
    }
    0
}

/// `UnMarkOneComponent`
pub fn unmark_one_component(at: &mut [InpAtom]) -> i32 {
    for a in at.iter_mut() {
        a.orig_compt_at_numb = 0;
        a.component = 0;
    }
    0
}

/// `UnMarkRingSystemsInp`
pub fn unmark_ring_systems_inp(at: &mut [InpAtom]) -> i32 {
    for a in at.iter_mut() {
        a.b_cut_vertex = 0;
        a.n_ring_system = 0;
        a.n_num_at_in_ring_system = 0;
        a.n_block_system = 0;
    }
    0
}

/// `MarkRingSystemsInp`: finds blocks (`n_block_system`, cut vertices) and
/// ring systems (`n_ring_system`, `n_num_at_in_ring_system`; a chain atom
/// between bridges forms a one-atom system) by two depth-first searches from
/// `start`. Returns the number of ring systems of the second search or a
/// negative error code.
pub fn mark_ring_systems_inp(at: &mut [InpAtom], start: usize) -> i32 {
    let num_atoms = at.len();
    if start >= num_atoms {
        return CT_OVERFLOW;
    }
    let mut stack_atom = vec![0u16; num_atoms];
    let mut ring_stack = vec![0u16; num_atoms];
    let mut dfs_number = vec![0u16; num_atoms];
    let mut low_number = vec![0u16; num_atoms];
    let mut neigh_numb = vec![0i8; num_atoms];

    // ---- cut vertices and blocks ----
    let mut num_ring_systems: i32 = 0;
    let mut dfs: u16 = 1;
    let mut top_stack: i32 = 0;
    let mut top_ring: i32 = 0;
    low_number[start] = dfs;
    dfs_number[start] = dfs;
    stack_atom[0] = start as u16;
    ring_stack[0] = start as u16;
    let mut num_start_children = 0;

    loop {
        // advance
        loop {
            let i = stack_atom[top_stack as usize] as usize;
            let j = neigh_numb[i];
            if at[i].valence <= j {
                neigh_numb[i] = 0;
                break;
            }
            neigh_numb[i] += 1;
            let u = at[i].neighbor[j as usize] as usize;
            if dfs_number[u] == 0 {
                top_stack += 1;
                stack_atom[top_stack as usize] = u as u16;
                top_ring += 1;
                ring_stack[top_ring as usize] = u as u16;
                dfs += 1;
                dfs_number[u] = dfs;
                low_number[u] = dfs;
                num_start_children += i32::from(i == start);
            } else if (top_stack == 0 || u != stack_atom[top_stack as usize - 1] as usize)
                && dfs_number[u] < dfs_number[i]
                && low_number[i] > dfs_number[u]
            {
                low_number[i] = dfs_number[u];
            }
        }
        // back up
        let i = stack_atom[top_stack as usize] as usize;
        if i != start {
            let u = stack_atom[top_stack as usize - 1] as usize;
            if low_number[i] >= dfs_number[u] {
                num_ring_systems += 1;
                at[u].n_block_system = num_ring_systems as u16;
                if u != start || num_start_children > 1 {
                    at[u].b_cut_vertex += 1;
                }
                while top_ring >= 0 {
                    let j = ring_stack[top_ring as usize] as usize;
                    top_ring -= 1;
                    at[j].n_block_system = num_ring_systems as u16;
                    if i == j {
                        break;
                    }
                }
            } else if low_number[u] > low_number[i] {
                low_number[u] = low_number[i];
            }
        }
        top_stack -= 1;
        if top_stack < 0 {
            break;
        }
    }

    // ---- ring systems ----
    num_ring_systems = 0;
    dfs = 1;
    top_stack = 0;
    top_ring = 0;
    dfs_number.iter_mut().for_each(|x| *x = 0);
    neigh_numb.iter_mut().for_each(|x| *x = 0);
    low_number[start] = dfs;
    dfs_number[start] = dfs;
    stack_atom[0] = start as u16;
    ring_stack[0] = start as u16;

    loop {
        loop {
            let i = stack_atom[top_stack as usize] as usize;
            let j = neigh_numb[i];
            if at[i].valence <= j {
                neigh_numb[i] = 0;
                break;
            }
            neigh_numb[i] += 1;
            let u = at[i].neighbor[j as usize] as usize;
            if dfs_number[u] == 0 {
                top_stack += 1;
                stack_atom[top_stack as usize] = u as u16;
                top_ring += 1;
                ring_stack[top_ring as usize] = u as u16;
                dfs += 1;
                dfs_number[u] = dfs;
                low_number[u] = dfs;
            } else if (top_stack == 0 || u != stack_atom[top_stack as usize - 1] as usize)
                && dfs_number[u] < dfs_number[i]
                && low_number[i] > dfs_number[u]
            {
                low_number[i] = dfs_number[u];
            }
        }
        let i = stack_atom[top_stack as usize] as usize;
        if dfs_number[i] == low_number[i] {
            num_ring_systems += 1;
            let mut num_at_in_ring_system: u16 = 0;
            let mut j = top_ring;
            while j >= 0 {
                num_at_in_ring_system = num_at_in_ring_system.wrapping_add(1);
                if i == ring_stack[j as usize] as usize {
                    break;
                }
                j -= 1;
            }
            while top_ring >= 0 {
                let j = ring_stack[top_ring as usize] as usize;
                top_ring -= 1;
                at[j].n_ring_system = num_ring_systems as u16;
                at[j].n_num_at_in_ring_system = num_at_in_ring_system;
                if i == j {
                    break;
                }
            }
        } else if top_stack > 0 {
            let j = stack_atom[top_stack as usize - 1] as usize;
            if low_number[j] > low_number[i] {
                low_number[j] = low_number[i];
            }
        }
        top_stack -= 1;
        if top_stack < 0 {
            break;
        }
    }
    num_ring_systems
}
