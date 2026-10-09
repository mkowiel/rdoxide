//! Ported from RDKit `MolChemicalFeatures/FeatureParser.{h,cpp}`: the feature
//! definition file format (`AtomType` and `DefineFeature … EndFeature`
//! blocks).

use std::collections::BTreeMap;
use std::fmt;
use std::sync::Arc;

use super::{MolChemicalFeatureDef, MolChemicalFeatureFactory};
use crate::RDKitError;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};

const BLANKS: &[char] = &[' ', '\t', '\r', '\n'];

/// `FeatureFileParseException`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureFileParseError {
    pub line_no: u32,
    pub line: String,
    pub message: String,
}

impl fmt::Display for FeatureFileParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Error parsing feature file at line {}:\n{}\n",
            self.line_no, self.message
        )
    }
}

/// What can go wrong when reading feature definitions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FeatureError {
    Parse(FeatureFileParseError),
    Rdkit(RDKitError),
}

impl From<RDKitError> for FeatureError {
    fn from(e: RDKitError) -> Self {
        FeatureError::Rdkit(e)
    }
}

impl fmt::Display for FeatureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FeatureError::Parse(e) => e.fmt(f),
            FeatureError::Rdkit(e) => e.fmt(f),
        }
    }
}

impl std::error::Error for FeatureError {}

type Result<T> = std::result::Result<T, FeatureError>;

fn parse_error<T>(line_no: u32, line: &str, message: impl Into<String>) -> Result<T> {
    Err(FeatureError::Parse(FeatureFileParseError {
        line_no,
        line: line.to_string(),
        message: message.into(),
    }))
}

/// The `std::getline` state machine of the text being parsed: `eof` is set
/// when a read reaches the end of the data, `fail` when it extracted nothing.
struct LineStream<'a> {
    rest: &'a str,
    eof: bool,
    fail: bool,
    line_no: u32,
}

impl<'a> LineStream<'a> {
    fn new(text: &'a str) -> Self {
        LineStream {
            rest: text,
            eof: false,
            fail: false,
            line_no: 0,
        }
    }

    fn exhausted(&self) -> bool {
        self.eof || self.fail
    }

    fn get_line(&mut self) -> &'a str {
        if self.rest.is_empty() {
            self.eof = true;
            self.fail = true;
            return "";
        }
        match self.rest.find('\n') {
            Some(i) => {
                let line = &self.rest[..i];
                self.rest = &self.rest[i + 1..];
                line
            }
            None => {
                self.eof = true;
                std::mem::take(&mut self.rest)
            }
        }
    }

    /// `getNextLine`: the next logical line. Blank and `#` comment lines are
    /// skipped and a trailing backslash continues the line. Like upstream,
    /// `line` is left alone when the stream is already exhausted.
    fn next_line(&mut self, line: &mut String) {
        if self.exhausted() {
            return;
        }
        line.clear();
        let mut continuation = false;
        while !self.exhausted() {
            let tmp = self.get_line();
            self.line_no += 1;
            if tmp.is_empty() || tmp.starts_with('#') {
                continue;
            }
            let stripped = tmp.trim_end_matches(BLANKS);
            if stripped.is_empty() {
                continue;
            }
            if !stripped.ends_with('\\') {
                line.push_str(if continuation {
                    tmp.trim_matches(BLANKS)
                } else {
                    tmp
                });
                return;
            }
            continuation = true;
            let trimmed = tmp.trim_matches(BLANKS);
            line.push_str(&trimmed[..trimmed.len() - 1]);
        }
    }
}

/// `boost::tokenizer` with `char_separator(" \t")`.
fn blank_tokens(line: &str) -> impl Iterator<Item = &str> {
    line.split([' ', '\t']).filter(|t| !t.is_empty())
}

/// The first token of boost's default tokenizer: a run of characters that are
/// neither blank nor punctuation, or a single punctuation character.
fn first_default_token(line: &str) -> Option<&str> {
    let line = line.trim_start_matches(|c: char| c.is_ascii_whitespace());
    let first = line.chars().next()?;
    if first.is_ascii_punctuation() {
        return Some(&line[..1]);
    }
    let end = line
        .find(|c: char| c.is_ascii_whitespace() || c.is_ascii_punctuation())
        .unwrap_or(line.len());
    Some(&line[..end])
}

/// `expandAndTestSmarts`: substitutes the atom-type labels and checks that
/// the result parses.
fn expand_and_test_smarts(smarts: &mut String, atom_types: &BTreeMap<String, String>) -> bool {
    for (name, definition) in atom_types {
        *smarts = smarts.replace(name.as_str(), definition);
    }
    matches!(
        mol_from_smarts(smarts, &SmartsParserParams::default()),
        Ok(Some(_))
    )
}

fn parse_atom_type(
    line: &str,
    atom_types: &mut BTreeMap<String, String>,
    line_no: u32,
) -> Result<()> {
    let mut tokens = blank_tokens(line);
    let Some(keyword) = tokens.next() else {
        return parse_error(line_no, line, "empty input line for AtomType");
    };
    if !keyword.eq_ignore_ascii_case("ATOMTYPE") {
        return parse_error(line_no, line, "bad input line for AtomType");
    }
    let Some(label) = tokens.next() else {
        return parse_error(line_no, line, "bad AtomType line, missing label");
    };
    let (negater, label) = match label.strip_prefix('!') {
        Some(rest) => (true, rest),
        None => (false, label),
    };
    let atom_type = format!("{{{label}}}");
    let Some(definition) = tokens.next() else {
        return parse_error(line_no, line, "bad AtomType line, missing definition");
    };
    let recursive = format!("$({definition})");
    let mut sma = match atom_types.get(&atom_type) {
        // extend the existing definition of this label
        Some(base) if negater => base.replacen('[', &format!("[!{recursive};"), 1),
        Some(base) => match base.rfind(']') {
            Some(i) => format!("{},{recursive}{}", &base[..i], &base[i..]),
            None => base.clone(),
        },
        None => recursive,
    };
    sma = format!("[{sma}]");
    if !expand_and_test_smarts(&mut sma, atom_types) {
        let msg = format!("invalid SMARTS in AtomType ({atom_type}): {sma}");
        return parse_error(line_no, line, msg);
    }
    // cut the brackets back off
    atom_types.insert(atom_type, sma[1..sma.len() - 1].to_string());
    Ok(())
}

fn parse_feature_def(
    stream: &mut LineStream<'_>,
    in_line: &str,
    atom_types: &BTreeMap<String, String>,
) -> Result<MolChemicalFeatureDef> {
    let mut next = in_line.trim_matches(BLANKS).to_string();
    while next.is_empty() || next.starts_with('#') {
        stream.next_line(&mut next);
        if next.is_empty() {
            return parse_error(
                stream.line_no,
                in_line,
                "EOF hit parsing feature definition",
            );
        }
        next = next.trim_matches(BLANKS).to_string();
    }
    let (sub_type, mut pattern) = {
        let mut tokens = blank_tokens(&next);
        if tokens.next().is_none() {
            return parse_error(
                stream.line_no,
                in_line,
                "bad DefineFeature line, no tokens found",
            );
        }
        let Some(sub_type) = tokens.next() else {
            return parse_error(
                stream.line_no,
                in_line,
                "bad DefineFeature line, missing subtype",
            );
        };
        let Some(pattern) = tokens.next() else {
            return parse_error(
                stream.line_no,
                in_line,
                "bad DefineFeature line, missing pattern",
            );
        };
        (sub_type.to_string(), pattern.to_string())
    };
    if !expand_and_test_smarts(&mut pattern, atom_types) {
        let msg = format!("invalid SMARTS in DefineFeature for type {sub_type}: {pattern}");
        return parse_error(stream.line_no, in_line, msg);
    }

    let mut weights = Vec::new();
    let mut family = String::new();
    let mut found_end = false;
    stream.next_line(&mut next);
    while !next.is_empty() {
        next = next.trim_matches(BLANKS).to_string();
        if !next.is_empty() && !next.starts_with('#') {
            let mut tokens = blank_tokens(&next);
            let keyword = tokens.next().unwrap_or_default().to_ascii_uppercase();
            match keyword.as_str() {
                "ENDFEATURE" => {
                    found_end = true;
                    break;
                }
                "FAMILY" => match tokens.next() {
                    Some(f) => family = f.to_string(),
                    None => {
                        let msg = format!("bad Type line for feature: {sub_type}");
                        return parse_error(stream.line_no, in_line, msg);
                    }
                },
                "WEIGHTS" => {
                    let Some(list) = tokens.next() else {
                        let msg = format!("bad Weights line for feature: {sub_type}");
                        return parse_error(stream.line_no, in_line, msg);
                    };
                    for number in list.split(',') {
                        match number.parse::<f64>() {
                            Ok(w) => weights.push(w),
                            Err(_) => {
                                let msg =
                                    format!("bad weight value ({number}) for feature: {sub_type}");
                                return parse_error(stream.line_no, in_line, msg);
                            }
                        }
                    }
                }
                _ => {
                    let msg = format!("bad input line for feature: {sub_type}");
                    return parse_error(stream.line_no, in_line, msg);
                }
            }
        }
        // upstream loops forever on a file that ends inside a definition;
        // report the missing EndFeature instead
        if stream.exhausted() {
            break;
        }
        stream.next_line(&mut next);
    }
    if !found_end {
        let msg = format!("could not find EndFeature line for feature: {sub_type}");
        return parse_error(stream.line_no, in_line, msg);
    }
    if family.is_empty() {
        let msg = format!("did not find Family definition for feature: {sub_type}");
        return parse_error(stream.line_no, in_line, msg);
    }

    let mut def = MolChemicalFeatureDef::new(&pattern, &family, &sub_type)?;
    if !weights.is_empty() {
        def.set_weights(&weights)?;
        def.normalize_weights();
    }
    Ok(def)
}

/// `parseFeatureData`: the feature definitions in `text`.
pub fn parse_feature_data(text: &str) -> Result<Vec<MolChemicalFeatureDef>> {
    let mut stream = LineStream::new(text);
    let mut defs = Vec::new();
    let mut atom_types = BTreeMap::new();
    let mut line = String::new();
    stream.next_line(&mut line);
    while !stream.exhausted() {
        line = line.trim_matches(BLANKS).to_string();
        if !line.is_empty() && !line.starts_with('#') {
            match first_default_token(&line)
                .map(str::to_ascii_uppercase)
                .as_deref()
            {
                Some("ATOMTYPE") => parse_atom_type(&line, &mut atom_types, stream.line_no)?,
                Some("DEFINEFEATURE") => {
                    let def = parse_feature_def(&mut stream, &line, &atom_types)?;
                    defs.push(def);
                }
                _ => return parse_error(stream.line_no, &line, "bad or missing keyword"),
            }
        }
        stream.next_line(&mut line);
    }
    Ok(defs)
}

/// `buildFeatureFactory`: a factory for the definitions in `text`.
pub fn build_feature_factory(text: &str) -> Result<MolChemicalFeatureFactory> {
    let mut factory = MolChemicalFeatureFactory::new();
    for def in parse_feature_data(text)? {
        factory.add_feature_def(Arc::new(def));
    }
    Ok(factory)
}

#[cfg(test)]
mod tests {
    use super::*;

    const FDEF: &str = "# comment\nAtomType NDonor [N&!H0&v3,N&!H0&+1&v4,n&H1&+0]\nAtomType AmideN [N;$(NC=O)]\nAtomType !NDonor [$(N(C)C)]\n\nDefineFeature HDonor1 [{NDonor}]\n  Family HBondDonor\n  Weights 1.0\nEndFeature\nDefineFeature Two [N;H2]\\\n   -C\n  Family X\n  Weights 1,3\nEndFeature\n";

    #[test]
    fn parses_atom_types_and_features() {
        let defs = parse_feature_data(FDEF).unwrap();
        assert_eq!(defs.len(), 2);
        assert_eq!(defs[0].family(), "HBondDonor");
        assert_eq!(defs[0].feature_type(), "HDonor1");
        assert_eq!(
            defs[0].smarts(),
            "[$([!$([$(N(C)C)]);N&!H0&v3,N&!H0&+1&v4,n&H1&+0])]"
        );
        assert_eq!(defs[1].smarts(), "[N;H2]-C");
        assert_eq!(defs[1].weights(), &[0.25, 0.75]);
    }

    #[test]
    fn errors_name_the_line() {
        for (text, needle) in [
            ("Bogus x y\n", "bad or missing keyword"),
            ("AtomType A\n", "missing definition"),
            ("AtomType A [N\n", "invalid SMARTS in AtomType"),
            (
                "DefineFeature a [N]\nFamily F\n",
                "could not find EndFeature",
            ),
            ("DefineFeature a [N]\nEndFeature\n", "did not find Family"),
            (
                "DefineFeature a [N]\nFamily F\nWeights x\nEndFeature\n",
                "bad weight value (x)",
            ),
            (
                "DefineFeature a [N]\nFamily F\nColor red\nEndFeature\n",
                "bad input line",
            ),
        ] {
            match parse_feature_data(text) {
                Err(FeatureError::Parse(e)) => assert!(e.message.contains(needle), "{text}: {e}"),
                other => panic!("{text}: {other:?}"),
            }
        }
        assert!(matches!(
            parse_feature_data("DefineFeature a [N]\nFamily F\nWeights 1,2\nEndFeature\n"),
            Err(FeatureError::Rdkit(RDKitError::Value(_)))
        ));
    }
}
