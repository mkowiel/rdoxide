//! Port of RDKit `Code/GraphMol/MolProcessing` and of the file-format
//! dispatch of `FileParsers/GeneralFileReader.h`: fingerprints for all the
//! molecules of a file.
//!
//! Not supported: TDT and Maestro files, compressed files (the RDKit builds
//! that lack Boost.IOStreams do not support them either).

use super::file_parsers::mol_file_parser::MolFileParserParams;
use super::file_parsers::multithreaded_supplier::{MultithreadedMolSupplier, SupplierParams};
use super::file_parsers::smiles_mol_supplier::SmilesMolSupplierParams;
use super::fingerprints::morgan::MorganArguments;
use super::fingerprints::{
    AtomInvariantsGenerator, BondInvariantsGenerator, EnvGenerator, FingerprintArguments,
    FingerprintFuncArguments, FingerprintGenerator,
};
use super::romol::RWMol;
use crate::datastructs::ExplicitBitVect;
use crate::rdgeneral::threads::par_map;
use crate::{RDKitError, Result};

/// `GeneralMolSupplier::SupplierOptions`
#[derive(Debug, Clone)]
pub struct SupplierOptions {
    pub sanitize: bool,
    pub remove_hs: bool,
    pub strict_parsing: bool,
    /// used for SMILES files
    pub delimiter: String,
    pub smiles_column: i32,
    pub name_column: i32,
    pub title_line: bool,
    /// used for TDT files
    pub name_record: String,
    pub conf_id_2d: i32,
    pub conf_id_3d: i32,
    /// the number of threads to use (<= 0: all cores minus `|n|`)
    pub num_writer_threads: i32,
}

impl Default for SupplierOptions {
    fn default() -> Self {
        SupplierOptions {
            sanitize: true,
            remove_hs: true,
            strict_parsing: true,
            delimiter: "\t".into(),
            smiles_column: 0,
            name_column: 1,
            title_line: true,
            name_record: String::new(),
            conf_id_2d: -1,
            conf_id_3d: 0,
            num_writer_threads: 0,
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
enum FileFormat {
    Sdf,
    Smiles,
    Tdt,
    Mae,
}

/// `determineFormat`: (format, compressed)
fn determine_format(path: &str) -> Result<(FileFormat, bool)> {
    let lower = path.to_lowercase();
    let unsupported = || {
        RDKitError::BadFile(format!(
            "Unsupported structure or compression extension given path: {path}"
        ))
    };
    if lower.ends_with(".maegz") {
        return Ok((FileFormat::Mae, true));
    }
    if lower.ends_with(".sdfgz") {
        return Ok((FileFormat::Sdf, true));
    }
    let (base, compressed) = if let Some(b) = lower.strip_suffix(".gz") {
        (b, true)
    } else if [".zst", ".bz2", ".7z"].iter().any(|e| lower.ends_with(e)) {
        return Err(RDKitError::BadFile(format!(
            "Unsupported compression extension (.zst, .bz2, .7z) given path: {path}"
        )));
    } else {
        (lower.as_str(), false)
    };
    let format = match base.rsplit_once('.').map(|(_, ext)| ext) {
        Some("sdf") => FileFormat::Sdf,
        Some("smi" | "csv" | "txt" | "tsv") => FileFormat::Smiles,
        Some("tdt") => FileFormat::Tdt,
        Some("mae") => FileFormat::Mae,
        _ => return Err(unsupported()),
    };
    Ok((format, compressed))
}

/// `GeneralMolSupplier::getSupplier`
fn get_supplier(path: &str, opt: &SupplierOptions) -> Result<MultithreadedMolSupplier> {
    let (format, compressed) = determine_format(path)?;
    if compressed {
        return Err(RDKitError::BadFile(
            "compressed files are only supported if the RDKit is built with boost::iostreams support"
                .into(),
        ));
    }
    let meta = std::fs::metadata(path)
        .map_err(|_| RDKitError::BadFile(format!("Bad input file {path}")))?;
    if meta.len() == 0 {
        return Err(RDKitError::BadFile(format!("Invalid input file {path}")));
    }
    let params = SupplierParams {
        num_writer_threads: opt.num_writer_threads,
        ..SupplierParams::default()
    };
    match format {
        FileFormat::Sdf => MultithreadedMolSupplier::from_sd_file(
            path,
            MolFileParserParams {
                sanitize: opt.sanitize,
                remove_hs: opt.remove_hs,
                strict_parsing: opt.strict_parsing,
                ..MolFileParserParams::default()
            },
            params,
        ),
        FileFormat::Smiles => {
            let mut parse_params = SmilesMolSupplierParams {
                delimiter: opt.delimiter.clone(),
                smiles_column: opt.smiles_column,
                name_column: opt.name_column,
                title_line: opt.title_line,
                ..SmilesMolSupplierParams::default()
            };
            parse_params.parse_parameters.sanitize = opt.sanitize;
            MultithreadedMolSupplier::from_smiles_file(path, parse_params, params)
        }
        FileFormat::Tdt | FileFormat::Mae => Err(RDKitError::BadFile(format!(
            "Unsupported file format: {}",
            if format == FileFormat::Tdt {
                "tdt"
            } else {
                "mae"
            }
        ))),
    }
}

/// `MorganFingerprint::getMorganGenerator(radius)` with its default arguments.
pub fn default_morgan_generator(radius: u32) -> FingerprintGenerator {
    FingerprintGenerator {
        env: EnvGenerator::Morgan(MorganArguments {
            radius,
            only_nonzero_invariants: false,
            include_redundant_environments: false,
        }),
        args: FingerprintArguments {
            count_simulation: false,
            include_chirality: false,
            count_bounds: vec![1, 2, 4, 8],
            fp_size: 2048,
            num_bits_per_feature: 1,
        },
        atom_inv_gen: Some(AtomInvariantsGenerator::Morgan {
            include_ring_membership: true,
        }),
        bond_inv_gen: Some(BondInvariantsGenerator::Morgan {
            use_bond_types: true,
            use_chirality: false,
        }),
        output_64: true,
    }
}

/// How many molecules are read before their fingerprints are computed.
const CHUNK: usize = 1024;

/// `getFingerprintsForMolsInFile`: the fingerprints of the molecules of a
/// file, in input order (`None` for the molecules that could not be read).
/// Without a `generator`, Morgan fingerprints of radius 3 are made.
pub fn get_fingerprints_for_mols_in_file(
    file_name: &str,
    options: &SupplierOptions,
    generator: Option<&FingerprintGenerator>,
) -> Result<Vec<Option<ExplicitBitVect>>> {
    let mut supplier = get_supplier(file_name, options)?;
    let default_generator;
    let generator = match generator {
        Some(g) => g,
        None => {
            default_generator = default_morgan_generator(3);
            &default_generator
        }
    };
    let mut results = Vec::new();
    loop {
        let mut chunk: Vec<Option<RWMol>> = Vec::with_capacity(CHUNK);
        while chunk.len() < CHUNK {
            match supplier.next_mol() {
                Some(mol) => chunk.push(mol),
                None => break,
            }
        }
        if chunk.is_empty() {
            break;
        }
        let last = chunk.len() < CHUNK;
        let fps = par_map(&chunk, options.num_writer_threads, |mol| {
            mol.as_ref()
                .map(|m| {
                    generator.fingerprint(
                        m,
                        &mut FingerprintFuncArguments {
                            from_atoms: None,
                            ignore_atoms: None,
                            conf_id: -1,
                            additional_output: None,
                            custom_atom_invariants: None,
                            custom_bond_invariants: None,
                        },
                    )
                })
                .transpose()
        });
        for fp in fps {
            results.push(fp?);
        }
        if last {
            break;
        }
    }
    // upstream reports no entries for unreadable molecules at the end of the
    // input (including the empty record that follows the last newline)
    while results.last().is_some_and(Option::is_none) {
        results.pop();
    }
    Ok(results)
}
