//! Ported from RDKit `Code/GraphMol/FilterCatalog/FunctionalGroupHierarchy.cpp`:
//! the built-in functional group definitions. Only the flattened name ->
//! pattern view is provided (the `FilterCatalog` form needs the catalog classes).

use std::collections::BTreeMap;
use std::sync::OnceLock;

use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::rdgeneral::RDProps;

/// One functional group: hierarchy `level`, dotted `name`, SMARTS, display
/// label and the optional reaction SMARTS removing the group.
#[allow(dead_code)] // `level` documents the hierarchy and is checked by the tests
struct FuncData {
    level: u8,
    name: &'static str,
    smarts: &'static str,
    label: &'static str,
    removal_reaction: Option<&'static str>,
}

const fn fd(
    level: u8,
    name: &'static str,
    smarts: &'static str,
    label: &'static str,
    removal_reaction: Option<&'static str>,
) -> FuncData {
    FuncData {
        level,
        name,
        smarts,
        label,
        removal_reaction,
    }
}

const FUNC_DATA: &[FuncData] = &[
    fd(0, "AcidChloride", "C(=O)Cl", "Acid Chloride", None),
    fd(
        1,
        "AcidChloride.Aromatic",
        "[$(C-!@[a])](=O)(Cl)",
        "Aromatic",
        None,
    ),
    fd(
        1,
        "AcidChloride.Aliphatic",
        "[$(C-!@[A;!Cl])](=O)(Cl)",
        "Aliphatic",
        None,
    ),
    fd(0, "CarboxylicAcid", "C(=O)[O;H,-]", "Carboxylic acid", None),
    fd(
        1,
        "CarboxylicAcid.Aromatic",
        "[$(C-!@[a])](=O)([O;H,-])",
        "Aromatic",
        None,
    ),
    fd(
        1,
        "CarboxylicAcid.Aliphatic",
        "[$(C-!@[A;!O])](=O)([O;H,-])",
        "Aliphatic",
        None,
    ),
    fd(
        1,
        "CarboxylicAcid.AlphaAmino",
        "[$(C-[C;!$(C=[!#6])]-[N;!H0;!$(N-[!#6;!#1]);!$(N-C=[O,N,S])])](=O)([O;H,-])",
        "alpha Amino Acid",
        None,
    ),
    fd(
        0,
        "SulfonylChloride",
        "[$(S-!@[#6])](=O)(=O)(Cl)",
        "Sulfonyl Chloride",
        None,
    ),
    fd(
        1,
        "SulfonylChloride.Aromatic",
        "[$(S-!@c)](=O)(=O)(Cl)",
        "Aromatic",
        None,
    ),
    fd(
        1,
        "SulfonylChloride.Aliphatic",
        "[$(S-!@C)](=O)(=O)(Cl)",
        "Aliphatic",
        None,
    ),
    fd(
        0,
        "Amine",
        "[N;$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]",
        "Amine",
        None,
    ),
    fd(
        1,
        "Amine.Primary",
        "[N;H2;D1;$(N-!@[#6]);!$(N-C=[O,N,S])]",
        "Primary",
        None,
    ),
    fd(
        2,
        "Amine.Primary.Aromatic",
        "[N;H2;D1;$(N-!@c);!$(N-C=[O,N,S])]",
        "Primary aromatic",
        None,
    ),
    fd(
        2,
        "Amine.Primary.Aliphatic",
        "[N;H2;D1;$(N-!@C);!$(N-C=[O,N,S])]",
        "Primary aliphatic",
        None,
    ),
    fd(
        1,
        "Amine.Secondary",
        "[N;H1;D2;$(N(-[#6])-[#6]);!$(N-C=[O,N,S])]",
        "Secondary",
        None,
    ),
    fd(
        2,
        "Amine.Secondary.Aromatic",
        "[N;H1;D2;$(N(-[c])-[#6]);!$(N-C=[O,N,S])]",
        "Secondary aromatic",
        None,
    ),
    fd(
        2,
        "Amine.Secondary.Aliphatic",
        "[N;H1;D2;$(N(-C)-C);!$(N-C=[O,N,S])]",
        "Secondary aliphatic",
        None,
    ),
    fd(
        1,
        "Amine.Tertiary",
        "[N;H0;D3;$(N(-[#6])(-[#6])-[#6]);!$(N-C=[O,N,S])]",
        "Tertiary",
        None,
    ),
    fd(
        2,
        "Amine.Tertiary.Aromatic",
        "[N;H0;D3;$(N(-[c])(-[#6])-[#6]);$(N-C=[O,N,S])]",
        "Tertiary aromatic",
        None,
    ),
    fd(
        2,
        "Amine.Tertiary.Aliphatic",
        "[N;H0;D3;$(N(-C)(-C)-C);!$(N-C=[O,N,S])]",
        "Tertiary aliphatic",
        None,
    ),
    fd(
        1,
        "Amine.Aromatic",
        "[N;$(N-c);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]",
        "Aromatic",
        None,
    ),
    fd(
        1,
        "Amine.Aliphatic",
        "[N;!$(N-c);$(N-C);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]",
        "Aliphatic",
        None,
    ),
    fd(
        1,
        "Amine.Cyclic",
        "[N;R;$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]",
        "Cyclic",
        None,
    ),
    fd(
        0,
        "BoronicAcid",
        "[$(B-!@[#6])](O)(O)",
        "Boronic Acid",
        Some("[#6:1]-!@[B:2]([O:3])[O:4]>>[#6:1][X].[B:2]([O:3])[O:4]"),
    ),
    fd(
        1,
        "BoronicAcid.Aromatic",
        "[$(B-!@c)](O)(O)",
        "Aromatic",
        Some("[c:1]-!@[B:2]([O:3])[O:4]>>[c:1][X].[B:2]([O:3])[O:4]"),
    ),
    fd(
        1,
        "BoronicAcid.Aliphatic",
        "[$(B-!@C)](O)(O)",
        "Aliphatic",
        Some("[C:1]-!@[B:2]([O:3])[O:4]>>[C:1][X].[B:2]([O:3])[O:4]"),
    ),
    fd(
        0,
        "Isocyanate",
        "[$(N-!@[#6])](=!@C=!@O)",
        "Isocyanate",
        None,
    ),
    fd(
        1,
        "Isocyanate.Aromatic",
        "[$(N-!@c)](=!@C=!@O)",
        "Aromatic",
        None,
    ),
    fd(
        1,
        "Isocyanate.Aliphatic",
        "[$(N-!@C)](=!@C=!@O)",
        "Aliphatic",
        None,
    ),
    fd(
        0,
        "Alcohol",
        "[O;H1;$(O-!@[#6;!$(C=!@[O,N,S])])]",
        "Alcohol",
        None,
    ),
    fd(1, "Alcohol.Aromatic", "[O;H1;$(O-!@c)]", "Aromatic", None),
    fd(
        1,
        "Alcohol.Aliphatic",
        "[O;H1;$(O-!@[C;!$(C=!@[O,N,S])])]",
        "Aliphatic",
        None,
    ),
    fd(0, "Aldehyde", "[CH;D2;!$(C-[!#6;!#1])]=O", "Aldehyde", None),
    fd(
        1,
        "Aldehyde.Aromatic",
        "[CH;D2;$(C-!@[a])](=O)",
        "Aromatic",
        None,
    ),
    fd(
        1,
        "Aldehyde.Aliphatic",
        "[CH;D2;$(C-!@C)](=O)",
        "Aliphatic",
        None,
    ),
    fd(
        0,
        "Halogen",
        "[$([F,Cl,Br,I]-!@[#6]);!$([F,Cl,Br,I]-!@C-!@[F,Cl,Br,I]);!$([F,Cl,Br,I]-[C,S](=[O,S,N]))]",
        "Halogen",
        None,
    ),
    fd(
        1,
        "Halogen.Aromatic",
        "[F,Cl,Br,I;$(*-!@c)]",
        "Aromatic",
        None,
    ),
    fd(
        1,
        "Halogen.Aliphatic",
        "[$([F,Cl,Br,I]-!@C);!$([F,Cl,Br,I]-!@C-!@[F,Cl,Br,I])]",
        "Aliphatic",
        None,
    ),
    fd(
        1,
        "Halogen.NotFluorine",
        "[$([Cl,Br,I]-!@[#6]);!$([Cl,Br,I]-!@C-!@[F,Cl,Br,I]);!$([Cl,Br,I]-[C,S](=[O,S,N]))]",
        "Not Fluorine",
        None,
    ),
    fd(
        2,
        "Halogen.NotFluorine.Aliphatic",
        "[$([Cl,Br,I]-!@C);!$([Cl,Br,I]-!@C-!@[F,Cl,Br,I]);!$([Cl,Br,I]-[C,S](=[O,S,N]))]",
        "Aliphatic Not Fluorine",
        None,
    ),
    fd(
        2,
        "Halogen.NotFluorine.Aromatic",
        "[$([Cl,Br,I]-!@c)]",
        "Aromatic Not Fluorine",
        None,
    ),
    fd(
        1,
        "Halogen.Bromine",
        "[Br;$([Br]-!@[#6]);!$([Br]-!@C-!@[F,Cl,Br,I]);!$([Br]-[C,S](=[O,S,N]))]",
        "Bromine",
        None,
    ),
    fd(
        2,
        "Halogen.Bromine.Aliphatic",
        "[Br;$(Br-!@C);!$(Br-!@C-!@[F,Cl,Br,I]);!$(Br-[C,S](=[O,S,N]))]",
        "Aliphatic Bromine",
        None,
    ),
    fd(
        2,
        "Halogen.Bromine.Aromatic",
        "[Br;$(Br-!@c)]",
        "Aromatic Bromine",
        None,
    ),
    fd(
        2,
        "Halogen.Bromine.BromoKetone",
        "[Br;$(Br-[CH2]-C(=O)-[#6])]",
        "Bromoketone",
        None,
    ),
    fd(
        0,
        "Azide",
        "[N;H0;$(N-[#6]);D2]=[N;D2]=[N;D1]",
        "Azide",
        None,
    ),
    fd(
        1,
        "Azide.Aromatic",
        "[N;H0;$(N-c);D2]=[N;D2]=[N;D1]",
        "Aromatic Azide",
        None,
    ),
    fd(
        1,
        "Azide.Aliphatic",
        "[N;H0;$(N-C);D2]=[N;D2]=[N;D1]",
        "Aliphatic Azide",
        None,
    ),
    fd(
        0,
        "Nitro",
        "[N;H0;$(N-[#6]);D3](=[O;D1])~[O;D1]",
        "Nitro",
        None,
    ),
    fd(
        1,
        "Nitro.Aromatic",
        "[N;H0;$(N-c);D3](=[O;D1])~[O;D1]",
        "Aromatic Nitro",
        None,
    ),
    fd(
        1,
        "Nitro.Aliphatic",
        "[N;H0;$(N-C);D3](=[O;D1])~[O;D1]",
        "Aliphatic Nitro",
        None,
    ),
    fd(
        0,
        "TerminalAlkyne",
        "[C;$(C#[CH])]",
        "Terminal Alkyne",
        None,
    ),
];

fn build(normalize: bool) -> BTreeMap<String, ROMol> {
    let params = SmartsParserParams {
        merge_hs: true,
        ..SmartsParserParams::default()
    };
    FUNC_DATA
        .iter()
        .map(|f| {
            let mut pattern = mol_from_smarts(f.smarts, &params)
                .ok()
                .flatten()
                .unwrap_or_else(|| panic!("bad functional group SMARTS for {}", f.name));
            pattern.dict_mut().set_val("Label", f.label);
            if let Some(reaction) = f.removal_reaction {
                pattern.dict_mut().set_val("RemovalReaction", reaction);
            }
            let key = if normalize {
                f.name.to_lowercase()
            } else {
                f.name.to_string()
            };
            (key, pattern)
        })
        .collect()
}

/// `GetFlattenedFunctionalGroupHierarchy(normalize)`: pattern molecules by
/// (optionally lower-cased) dotted name.
pub fn flattened_functional_group_hierarchy(normalize: bool) -> &'static BTreeMap<String, ROMol> {
    static EXACT: OnceLock<BTreeMap<String, ROMol>> = OnceLock::new();
    static NORMALIZED: OnceLock<BTreeMap<String, ROMol>> = OnceLock::new();
    if normalize {
        NORMALIZED.get_or_init(|| build(true))
    } else {
        EXACT.get_or_init(|| build(false))
    }
}

/// The hierarchy depth of every entry is at most this deep.
pub const MAX_DEPTH: usize = 3;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_patterns_parse() {
        let flat = flattened_functional_group_hierarchy(false);
        assert_eq!(flat.len(), FUNC_DATA.len());
        assert!(flat.contains_key("CarboxylicAcid.Aromatic"));
        assert!(flattened_functional_group_hierarchy(true).contains_key("acidchloride"));
    }

    #[test]
    fn levels_match_dotted_names() {
        for f in FUNC_DATA {
            assert_eq!(
                f.name.matches('.').count(),
                usize::from(f.level),
                "{}",
                f.name
            );
            assert!((f.level as usize) < MAX_DEPTH);
        }
    }
}
