//! `cargo run --example draw_png -- SMILES out.png`: draws a molecule as PNG (a check of the PNG drawer)
use rdoxide::graphmol::depictor::{Compute2DCoordParameters, compute_2d_coords};
use rdoxide::graphmol::mol_draw2d::raster::MolDraw2DCairo;
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut mol = mol_from_smiles(&args[1], &SmilesParserParams::default())
        .unwrap()
        .unwrap();
    compute_2d_coords(&mut mol, &Compute2DCoordParameters::default()).unwrap();
    let mut d = MolDraw2DCairo::new(400, 300, -1, -1).unwrap();
    d.draw_molecule(&mol, "legend", -1).unwrap();
    d.finish_drawing();
    d.write_drawing_text(&args[2]).unwrap();
}
