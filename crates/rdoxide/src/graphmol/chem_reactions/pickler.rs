//! Ported from RDKit `Code/GraphMol/ChemReactions/ReactionPickler.{h,cpp}`:
//! the binary reaction format (templates are stored as molecule pickles).

use super::ChemicalReaction;
use crate::graphmol::pickle::flags::{ALL_PROPS, ATOM_PROPS, MOL_PROPS};
use crate::graphmol::pickle::{
    ENDIAN_ID, get_default_pickle_properties, mol_from_pickle_stream, pickle_mol,
    pickle_properties, unpickle_properties,
};
use crate::graphmol::romol::ROMol;
use crate::{RDKitError, Result, rd_log_warning};

const VERSION_MAJOR: i32 = 3;
const VERSION_MINOR: i32 = 0;
const VERSION_PATCH: i32 = 0;

/// `ReactionPickler::Tags`
mod tag {
    pub const VERSION: i32 = 10000;
    pub const BEGIN_REACTANTS: i32 = 10001;
    pub const END_REACTANTS: i32 = 10002;
    pub const BEGIN_PRODUCTS: i32 = 10003;
    pub const END_PRODUCTS: i32 = 10004;
    pub const BEGIN_AGENTS: i32 = 10005;
    pub const END_AGENTS: i32 = 10006;
    pub const END_REACTION: i32 = 10007;
    pub const BEGIN_PROPS: i32 = 10008;
    pub const END_PROPS: i32 = 10009;
    pub const BEGIN_SSS_PARAMS: i32 = 10010;
    pub const END_SSS_PARAMS: i32 = 10011;
}

/// `ReactionPicklerException` (a `RuntimeError` in Python)
fn pickle_error(msg: impl Into<String>) -> RDKitError {
    RDKitError::Runtime(msg.into())
}

fn write_i32(out: &mut Vec<u8>, v: i32) {
    out.extend_from_slice(&v.to_le_bytes());
}

fn write_u32(out: &mut Vec<u8>, v: u32) {
    out.extend_from_slice(&v.to_le_bytes());
}

/// Sequential little-endian reader over a pickle.
struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl Reader<'_> {
    fn u32(&mut self) -> Result<u32> {
        let end = self.pos + 4;
        let bytes = self.data.get(self.pos..end).ok_or_else(|| {
            pickle_error("Bad pickle format: unexpected End-of-File while reading")
        })?;
        self.pos = end;
        Ok(u32::from_le_bytes(bytes.try_into().expect("four bytes")))
    }

    fn i32(&mut self) -> Result<i32> {
        Ok(self.u32()? as i32)
    }

    /// Reads a tag and checks that it is `expected`.
    fn expect_tag(&mut self, expected: i32, name: &str) -> Result<()> {
        if self.i32()? == expected {
            Ok(())
        } else {
            Err(pickle_error(format!(
                "Bad pickle format: {name} tag not found."
            )))
        }
    }

    fn mol(&mut self) -> Result<ROMol> {
        mol_from_pickle_stream(self.data, &mut self.pos, ALL_PROPS)
    }

    fn mols(&mut self, n: u32, begin: (i32, &str), end: (i32, &str)) -> Result<Vec<ROMol>> {
        self.expect_tag(begin.0, begin.1)?;
        let mols = (0..n).map(|_| self.mol()).collect::<Result<Vec<_>>>()?;
        self.expect_tag(end.0, end.1)?;
        Ok(mols)
    }
}

/// `ReactionPickler::pickleReaction(rxn, res, propertyFlags)`
pub fn pickle_reaction(rxn: &ChemicalReaction, property_flags: u32) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    for v in [
        ENDIAN_ID,
        tag::VERSION,
        VERSION_MAJOR,
        VERSION_MINOR,
        VERSION_PATCH,
    ] {
        write_i32(&mut out, v);
    }
    write_u32(&mut out, rxn.num_reactant_templates() as u32);
    write_u32(&mut out, rxn.num_product_templates() as u32);
    write_u32(&mut out, rxn.num_agent_templates() as u32);
    write_u32(
        &mut out,
        u32::from(rxn.implicit_properties) | (u32::from(rxn.needs_init) << 1),
    );

    // an initialized reaction carries private properties set by `validate`
    let template_flags = if rxn.needs_init {
        property_flags
    } else {
        ALL_PROPS
    };
    write_i32(&mut out, tag::BEGIN_REACTANTS);
    for mol in &rxn.reactants {
        out.extend(pickle_mol(mol, template_flags)?);
    }
    write_i32(&mut out, tag::END_REACTANTS);
    write_i32(&mut out, tag::BEGIN_PRODUCTS);
    for mol in &rxn.products {
        out.extend(pickle_mol(mol, template_flags)?);
    }
    write_i32(&mut out, tag::END_PRODUCTS);
    if !rxn.agents.is_empty() {
        write_i32(&mut out, tag::BEGIN_AGENTS);
        for mol in &rxn.agents {
            out.extend(pickle_mol(mol, property_flags | ATOM_PROPS)?);
        }
        write_i32(&mut out, tag::END_AGENTS);
    }

    write_i32(&mut out, tag::BEGIN_SSS_PARAMS);
    let json = rxn.substruct_params.to_json();
    write_u32(&mut out, json.len() as u32);
    out.extend_from_slice(json.as_bytes());
    write_i32(&mut out, tag::END_SSS_PARAMS);

    if property_flags & MOL_PROPS != 0 {
        write_i32(&mut out, tag::BEGIN_PROPS);
        out.extend(pickle_properties(rxn, property_flags));
        write_i32(&mut out, tag::END_PROPS);
    }
    write_i32(&mut out, tag::END_REACTION);
    Ok(out)
}

/// `ReactionPickler::pickleReaction(rxn, res)` with the default property flags.
pub fn pickle_reaction_default(rxn: &ChemicalReaction) -> Result<Vec<u8>> {
    pickle_reaction(rxn, get_default_pickle_properties())
}

/// `ReactionPickler::reactionFromPickle`
pub fn reaction_from_pickle(data: &[u8]) -> Result<ChemicalReaction> {
    let mut r = Reader { data, pos: 0 };
    if r.i32()? != ENDIAN_ID {
        return Err(pickle_error(
            "Bad pickle format: bad endian ID or invalid file format",
        ));
    }
    r.expect_tag(tag::VERSION, "version")
        .map_err(|_| pickle_error("Bad pickle format: no version tag"))?;
    let (major, minor, patch) = (r.i32()?, r.i32()?, r.i32()?);
    if major > VERSION_MAJOR || (major == VERSION_MAJOR && minor > VERSION_MINOR) {
        rd_log_warning!(
            "Depickling from a version number ({major}.{minor})that is higher than our version ({VERSION_MAJOR}.{VERSION_MINOR}).\nThis probably won't work."
        );
    }
    let version = 1000 * major + minor * 10 + patch;

    let n_reactants = r.u32()?;
    let n_products = r.u32()?;
    let n_agents = if version > 1000 { r.u32()? } else { 0 };
    let flag = r.u32()?;

    let mut rxn = ChemicalReaction::new();
    rxn.set_implicit_properties_flag(flag & 0x1 != 0);
    for mol in r.mols(
        n_reactants,
        (tag::BEGIN_REACTANTS, "BEGINREACTANTS"),
        (tag::END_REACTANTS, "ENDREACTANTS"),
    )? {
        rxn.add_reactant_template(mol);
    }
    for mol in r.mols(
        n_products,
        (tag::BEGIN_PRODUCTS, "BEGINPRODUCTS"),
        (tag::END_PRODUCTS, "ENDPRODUCTS"),
    )? {
        rxn.add_product_template(mol);
    }
    if n_agents != 0 {
        for mol in r.mols(
            n_agents,
            (tag::BEGIN_AGENTS, "BEGINAGENTS"),
            (tag::END_AGENTS, "ENDAGENTS"),
        )? {
            rxn.add_agent_template(mol);
        }
    }

    let mut next = r.i32();
    if version >= 3000 && next.as_ref().is_ok_and(|&t| t == tag::BEGIN_SSS_PARAMS) {
        let len = r.u32()? as usize;
        let bytes = r.data.get(r.pos..r.pos + len).ok_or_else(|| {
            pickle_error("Bad pickle format: unexpected End-of-File while reading")
        })?;
        r.pos += len;
        rxn.substruct_params
            .update_from_json(&String::from_utf8_lossy(bytes))?;
        r.expect_tag(tag::END_SSS_PARAMS, "ENDSSSPARAMS")?;
        next = r.i32();
    }
    if next.as_ref().is_ok_and(|&t| t == tag::BEGIN_PROPS) {
        unpickle_properties(data, &mut r.pos, &mut rxn.props)?;
        r.expect_tag(tag::END_PROPS, "ENDPROPS")?;
    }
    // this has to be set after the templates were added
    rxn.needs_init = flag & 0x2 != 0;
    Ok(rxn)
}
