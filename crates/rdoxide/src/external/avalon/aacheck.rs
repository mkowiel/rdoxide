//! Port of Avalon Toolkit's `aacheck.c`: checks that every atom matches one of the allowed
//! augmented atoms.

use super::pattern::{AugmentedAtom, aa_match, aa_to_string};
use super::perceive::{CHAIN, RING, ring_state, setup_neighbourhood};
use super::reaccs::Molecule;
use super::stereo::add_msg;

/// `CheckAtoms`: true if every atom of `mp` matches one of `good_atoms`; messages name the
/// atoms that do not
pub fn check_atoms(
    mp: &Molecule,
    good_atoms: &[AugmentedAtom],
    msgs: &mut Vec<String>,
) -> Result<bool, &'static str> {
    let (atom_status, _) = ring_state(mp);
    let mut nbp = setup_neighbourhood(mp).ok_or("too many neighbours")?;
    let mut nmatch = 0;
    for i in 0..mp.atoms.len() {
        let mut found = false;
        for g in good_atoms {
            if g.topography == RING && atom_status[i] == 0 {
                continue;
            }
            if g.topography == CHAIN && atom_status[i] != 0 {
                continue;
            }
            if nbp[i].atoms.len() == g.ligands.len()
                && aa_match(mp, i, g, Some(&atom_status), &nbp).is_some()
            {
                nmatch += 1;
                found = true;
                break;
            }
        }
        if !found {
            let s = aa_to_string(mp, i, &mut nbp[i]);
            add_msg(
                msgs,
                format!("{:>10}    atom {:3}    AA : {}", mp.name, i + 1, s),
            );
        }
    }
    Ok(nmatch == mp.atoms.len())
}
