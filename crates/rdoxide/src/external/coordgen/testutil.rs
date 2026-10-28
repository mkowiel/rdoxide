//! helpers of the C++ oracle tests: molecules from the input blocks of `tools/coordgen_drv/*.cpp`

use super::maths::PointF;
use super::model::{Model, MolId};

/// "index list" text of the drivers: positions of `ids` in `all`
pub fn list(ids: &[usize], all: &[usize]) -> String {
    ids.iter()
        .map(|i| format!(" {}", all.iter().position(|x| x == i).unwrap()))
        .collect()
}

/// builds the molecule of a "MOL n m ... END" block (atoms may carry the flags `f`ixed and `c`onstrained)
pub fn build(input: &str) -> (Model, MolId) {
    let mut model = Model::new();
    let mol = model.add_molecule();
    for line in input.lines().skip(1) {
        let f: Vec<&str> = line.split(' ').collect();
        match f[0] {
            "A" => {
                let a = model.add_atom(mol);
                model.atoms[a].atomic_number = f[1].parse().unwrap();
                model.atoms[a].charge = f[2].parse().unwrap();
                let x = f32::from_bits(u32::from_str_radix(f[3], 16).unwrap());
                let y = f32::from_bits(u32::from_str_radix(f[4], 16).unwrap());
                model.set_coordinates(a, PointF::new(x, y));
                if let Some(flags) = f.get(5) {
                    model.atoms[a].fixed = flags.contains('f');
                    model.atoms[a].constrained = flags.contains('c');
                }
            }
            "B" => {
                let b = model.add_bond(mol, f[1].parse().unwrap(), f[2].parse().unwrap());
                model.bonds[b].bond_order = f[3].parse().unwrap();
            }
            _ => {}
        }
    }
    (model, mol)
}
