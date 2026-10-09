//! `rdkit.Chem.rdqueries` — port of upstream `Code/GraphMol/Wrap/queries.cpp`:
//! factories for `QueryAtom`/`QueryBond` objects and `ReplaceAtomWithQueryAtom`.
//!
//! The upstream module is a wall of macro-generated functions; here the
//! families are tables and each function is a closure registered under its
//! generated name.

use pyo3::exceptions::PyOverflowError;
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyCFunction, PyDict, PyFloat, PyInt, PyString, PyTuple};
use rdoxide::graphmol::queries::{
    AtomFn, AtomQuery, BondFn, BondQuery, MASS_INTEGER_CONVERSION_FACTOR, PropValue, Query,
    QueryKind, make_a_atom_query, make_ah_atom_query, make_atom_aliphatic_query,
    make_atom_aromatic_query, make_atom_explicit_degree_query, make_atom_explicit_valence_query,
    make_atom_formal_charge_query, make_atom_h_count_query, make_atom_has_chiral_tag_query,
    make_atom_hybridization_query, make_atom_in_n_rings_query, make_atom_in_ring_query,
    make_atom_is_bridgehead_query, make_atom_isotope_query, make_atom_mass_query,
    make_atom_min_ring_size_query, make_atom_missing_chiral_tag_query,
    make_atom_non_hydrogen_degree_query, make_atom_num_aliphatic_heteroatom_nbrs_query,
    make_atom_num_heteroatom_nbrs_query, make_atom_num_query,
    make_atom_num_radical_electrons_query, make_atom_ring_bond_count_query,
    make_atom_total_degree_query, make_atom_total_valence_query, make_atom_unsaturated_query,
    make_has_prop_query, make_m_atom_query, make_mh_atom_query, make_prop_query, make_q_atom_query,
    make_qh_atom_query, make_x_atom_query, make_xh_atom_query,
};
use rdoxide::graphmol::query_ops::replace_atom_with_query_atom;
use rdoxide::graphmol::{Atom, Bond};

use crate::atom::PyAtom;
use crate::bond::PyBond;
use crate::errors::{argument_error, boost_type_name};
use crate::mol::PyMol;
use crate::rdmolfiles_write::Args;

const STRING_TYPE: &str =
    "std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>>";

/// How a Boost.Python argument converter treats a Python object.
#[derive(Clone, Copy)]
enum Kind {
    Int,
    Bool,
    Double,
    Str,
}

impl Kind {
    fn cpp(self) -> &'static str {
        match self {
            Kind::Int => "int",
            Kind::Bool => "bool",
            Kind::Double => "double",
            Kind::Str => STRING_TYPE,
        }
    }
}

/// A converted argument.
#[derive(Clone)]
enum Value {
    Int(i32),
    Bool(bool),
    Double(f64),
    Str(String),
}

/// One parameter of a factory: its name, converter and default (if any).
struct Param {
    name: &'static str,
    kind: Kind,
    default: Option<Value>,
}

impl Param {
    fn required(name: &'static str, kind: Kind) -> Self {
        Param {
            name,
            kind,
            default: None,
        }
    }

    fn optional(name: &'static str, kind: Kind, default: Value) -> Self {
        Param {
            name,
            kind,
            default: Some(default),
        }
    }

    fn signature(&self) -> String {
        let default = match &self.default {
            None => String::new(),
            Some(Value::Bool(b)) => format!("={}", if *b { "True" } else { "False" }),
            Some(Value::Int(i)) => format!("={i}"),
            Some(Value::Double(d)) => format!("={d:?}"),
            Some(Value::Str(s)) => format!("={s}"),
        };
        format!("{} {}{default}", self.kind.cpp(), self.name)
    }

    /// `None` when the object is not convertible (an `ArgumentError` upstream).
    fn convert(&self, obj: &Bound<'_, PyAny>) -> PyResult<Option<Value>> {
        let is_int = obj.is_instance_of::<PyInt>();
        Ok(match self.kind {
            Kind::Int if is_int => match obj.extract::<i32>() {
                Ok(v) => Some(Value::Int(v)),
                Err(e) if e.is_instance_of::<PyOverflowError>(obj.py()) => return Err(e),
                Err(_) => None,
            },
            Kind::Bool if obj.is_none() => Some(Value::Bool(false)),
            Kind::Bool if is_int || obj.is_instance_of::<PyBool>() => {
                Some(Value::Bool(obj.is_truthy()?))
            }
            Kind::Double if is_int || obj.is_instance_of::<PyFloat>() => {
                Some(Value::Double(obj.extract()?))
            }
            Kind::Str if obj.is_instance_of::<PyString>() => Some(Value::Str(obj.extract()?)),
            _ => None,
        })
    }
}

/// The converted arguments of a call, in parameter order.
struct Parsed(Vec<Value>);

impl Parsed {
    fn int(&self, i: usize) -> i32 {
        match self.0[i] {
            Value::Int(v) => v,
            _ => unreachable!("parameter {i} is not an int"),
        }
    }
    fn bool(&self, i: usize) -> bool {
        matches!(self.0[i], Value::Bool(true))
    }
    fn double(&self, i: usize) -> f64 {
        match self.0[i] {
            Value::Double(v) => v,
            _ => unreachable!("parameter {i} is not a double"),
        }
    }
    fn str(&self, i: usize) -> String {
        match &self.0[i] {
            Value::Str(v) => v.clone(),
            _ => unreachable!("parameter {i} is not a string"),
        }
    }
}

fn leak(s: String) -> &'static str {
    Box::leak(s.into_boxed_str())
}

/// Registers `name(params...)` on `m`; `build` turns the converted arguments
/// into the returned object. Calls that do not fit raise Boost's
/// `ArgumentError`.
fn def_factory<F>(
    m: &Bound<'_, PyModule>,
    name: String,
    doc: String,
    params: Vec<Param>,
    build: F,
) -> PyResult<()>
where
    F: Fn(Python<'_>, &Parsed) -> PyResult<Py<PyAny>> + Send + Sync + 'static,
{
    let py = m.py();
    let static_name: &'static str = leak(name);
    let c_name: &'static std::ffi::CStr = leak_cstr(static_name);
    let c_doc: &'static std::ffi::CStr = leak_cstr(&doc);
    let signature = format!(
        "{static_name}({})",
        params
            .iter()
            .map(Param::signature)
            .collect::<Vec<_>>()
            .join(", ")
    );
    let closure = move |args: &Bound<'_, PyTuple>,
                        kwargs: Option<&Bound<'_, PyDict>>|
          -> PyResult<Py<PyAny>> {
        let py = args.py();
        let bad_call = || {
            let got: Vec<String> = args.iter().map(|a| boost_type_name(&a)).collect();
            argument_error(py, static_name, &got, &signature)
        };
        // a keyword must name a parameter that was not given positionally
        let unusable = |key: &str| {
            params
                .iter()
                .position(|p| p.name == key)
                .is_none_or(|i| i < args.len())
        };
        if args.len() > params.len()
            || kwargs.is_some_and(|k| {
                k.keys()
                    .iter()
                    .any(|key| key.extract::<String>().map_or(true, |key| unusable(&key)))
            })
        {
            return Err(bad_call());
        }
        let call = Args {
            name: static_name,
            args,
            kwargs,
        };
        let mut values = Vec::with_capacity(params.len());
        for (i, param) in params.iter().enumerate() {
            let value = match call.get(i, param.name)? {
                Some(obj) => param.convert(&obj)?,
                None => param.default.clone(),
            };
            values.push(value.ok_or_else(bad_call)?);
        }
        build(py, &Parsed(values))
    };
    m.add_function(PyCFunction::new_closure(
        py,
        Some(c_name),
        Some(c_doc),
        closure,
    )?)
}

fn leak_cstr(s: &str) -> &'static std::ffi::CStr {
    Box::leak(
        std::ffi::CString::new(s)
            .expect("no NUL in names")
            .into_boxed_c_str(),
    )
}

fn negate_param() -> Param {
    Param::optional("negate", Kind::Bool, Value::Bool(false))
}

fn atom_object(py: Python<'_>, mut query: AtomQuery, negate: bool) -> PyResult<Py<PyAny>> {
    // like upstream this only ever switches negation on
    query.negate |= negate;
    let mut atom = Atom::new(0);
    atom.set_query(Some(query));
    Ok(PyAtom::owned(py, atom)?.into_any())
}

fn bond_object(py: Python<'_>, mut query: BondQuery, negate: bool) -> PyResult<Py<PyAny>> {
    query.negate |= negate;
    let mut bond = Bond::new(rdoxide::graphmol::BondType::Unspecified);
    bond.set_query(Some(query));
    Ok(PyBond::owned(py, bond)?.into_any())
}

/// An `...EqualsQueryAtom`/`...LessQueryAtom`/`...GreaterQueryAtom` family.
/// The C++ comparison queries test `value <op> atom property`; upstream
/// deliberately swaps Less and Greater so that `AtomNumLessQueryAtom(7)`
/// matches atoms whose atomic number is less than 7.
struct Comparison {
    name: &'static str,
    /// the equality query, which has the family's specific shape and description
    equals: fn(i32) -> AtomQuery,
    func: AtomFn,
    /// factor applied to the value (the mass is stored as an integer)
    scale: i32,
}

const COMPARISONS: [Comparison; 17] = [
    Comparison {
        name: "AtomNum",
        equals: make_atom_num_query,
        func: AtomFn::Num,
        scale: 1,
    },
    Comparison {
        name: "ExplicitValence",
        equals: make_atom_explicit_valence_query,
        func: AtomFn::ExplicitValence,
        scale: 1,
    },
    Comparison {
        name: "TotalValence",
        equals: make_atom_total_valence_query,
        func: AtomFn::TotalValence,
        scale: 1,
    },
    Comparison {
        name: "ExplicitDegree",
        equals: make_atom_explicit_degree_query,
        func: AtomFn::ExplicitDegree,
        scale: 1,
    },
    Comparison {
        name: "TotalDegree",
        equals: make_atom_total_degree_query,
        func: AtomFn::TotalDegree,
        scale: 1,
    },
    Comparison {
        name: "HCount",
        equals: make_atom_h_count_query,
        func: AtomFn::HCount,
        scale: 1,
    },
    Comparison {
        name: "Mass",
        equals: make_atom_mass_query,
        func: AtomFn::Mass,
        scale: MASS_INTEGER_CONVERSION_FACTOR,
    },
    Comparison {
        name: "Isotope",
        equals: make_atom_isotope_query,
        func: AtomFn::Isotope,
        scale: 1,
    },
    Comparison {
        name: "FormalCharge",
        equals: make_atom_formal_charge_query,
        func: AtomFn::FormalCharge,
        scale: 1,
    },
    Comparison {
        name: "Hybridization",
        equals: make_atom_hybridization_query,
        func: AtomFn::Hybridization,
        scale: 1,
    },
    Comparison {
        name: "InNRings",
        equals: make_atom_in_n_rings_query,
        func: AtomFn::InNRings,
        scale: 1,
    },
    Comparison {
        name: "MinRingSize",
        equals: make_atom_min_ring_size_query,
        func: AtomFn::MinRingSize,
        scale: 1,
    },
    Comparison {
        name: "RingBondCount",
        equals: make_atom_ring_bond_count_query,
        func: AtomFn::RingBondCount,
        scale: 1,
    },
    Comparison {
        name: "NumRadicalElectrons",
        equals: make_atom_num_radical_electrons_query,
        func: AtomFn::NumRadicalElectrons,
        scale: 1,
    },
    Comparison {
        name: "NumHeteroatomNeighbors",
        equals: make_atom_num_heteroatom_nbrs_query,
        func: AtomFn::NumHeteroatomNbrs,
        scale: 1,
    },
    Comparison {
        name: "NumAliphaticHeteroatomNeighbors",
        equals: make_atom_num_aliphatic_heteroatom_nbrs_query,
        func: AtomFn::NumAliphaticHeteroatomNbrs,
        scale: 1,
    },
    Comparison {
        name: "NonHydrogenDegree",
        equals: make_atom_non_hydrogen_degree_query,
        func: AtomFn::NonHydrogenDegree,
        scale: 1,
    },
];

/// Queries that are true or false for an atom.
type FlagQuery = (&'static str, fn() -> AtomQuery);

const FLAGS: [FlagQuery; 15] = [
    ("IsUnsaturated", make_atom_unsaturated_query),
    ("IsAromatic", make_atom_aromatic_query),
    ("IsAliphatic", make_atom_aliphatic_query),
    ("IsInRing", make_atom_in_ring_query),
    ("HasChiralTag", make_atom_has_chiral_tag_query),
    ("MissingChiralTag", make_atom_missing_chiral_tag_query),
    ("IsBridgehead", make_atom_is_bridgehead_query),
    ("AAtom", make_a_atom_query),
    ("AHAtom", make_ah_atom_query),
    ("XAtom", make_x_atom_query),
    ("XHAtom", make_xh_atom_query),
    ("QAtom", make_q_atom_query),
    ("QHAtom", make_qh_atom_query),
    ("MAtom", make_m_atom_query),
    ("MHAtom", make_mh_atom_query),
];

fn def_comparisons(m: &Bound<'_, PyModule>) -> PyResult<()> {
    for cmp in COMPARISONS {
        let params = || vec![Param::required("val", Kind::Int), negate_param()];
        let Comparison {
            name,
            equals,
            func,
            scale,
        } = cmp;
        def_factory(
            m,
            format!("{name}EqualsQueryAtom"),
            format!(
                "Returns a QueryAtom that matches atoms where {name} is equal to the target value."
            ),
            params(),
            move |py, a| atom_object(py, equals(a.int(0)), a.bool(1)),
        )?;
        // the direction of comparison is reversed relative to the C++ API
        for (suffix, kind, doc) in [
            ("Less", QueryKind::Greater, "is less than the target value."),
            ("Greater", QueryKind::Less, "is equal to the target value."),
        ] {
            def_factory(
                m,
                format!("{name}{suffix}QueryAtom"),
                format!(
                    "Returns a QueryAtom that matches atoms where {name} {doc}\n\
                     NOTE: the direction of comparison is reversed relative to the C++ API"
                ),
                params(),
                {
                    let func = func.clone();
                    move |py: Python<'_>, a: &Parsed| {
                        let query = Query::with(
                            kind.clone(),
                            scale * a.int(0),
                            func.clone(),
                            &format!("{name}{suffix}"),
                        );
                        atom_object(py, query, a.bool(1))
                    }
                },
            )?;
        }
    }
    Ok(())
}

fn def_flags(m: &Bound<'_, PyModule>) -> PyResult<()> {
    for (name, make) in FLAGS {
        let label = name.strip_suffix("Atom").unwrap_or(name);
        def_factory(
            m,
            format!("{name}QueryAtom"),
            format!("Returns a QueryAtom that matches atoms when {label} is True."),
            vec![negate_param()],
            move |py, a| atom_object(py, make(), a.bool(0)),
        )?;
    }
    Ok(())
}

/// The `Has...PropWithValueQuery{Atom,Bond}` families.
fn def_property_queries(m: &Bound<'_, PyModule>) -> PyResult<()> {
    def_factory(
        m,
        "HasPropQueryAtom".into(),
        "Returns a QueryAtom that matches when the property 'propname' exists in the atom.".into(),
        vec![Param::required("propname", Kind::Str), negate_param()],
        |py, a| atom_object(py, make_has_prop_query(&a.str(0)), a.bool(1)),
    )?;
    def_factory(
        m,
        "HasPropQueryBond".into(),
        "Returns a QueryBond that matches when the property 'propname' exists in the bond.".into(),
        vec![Param::required("propname", Kind::Str), negate_param()],
        |py, a| bond_object(py, make_has_prop_query(&a.str(0)), a.bool(1)),
    )?;
    type Value = fn(&Parsed) -> (PropValue, f64);
    let kinds: [(&str, &str, Kind, Option<Param>, Value); 4] = [
        (
            "Int",
            "int value",
            Kind::Int,
            Some(Param::optional("tolerance", Kind::Int, self::Value::Int(0))),
            |a| (PropValue::Int(a.int(1)), f64::from(a.int(3))),
        ),
        ("Bool", "boolean value", Kind::Bool, None, |a| {
            (PropValue::Bool(a.bool(1)), 0.0)
        }),
        ("String", "string value", Kind::Str, None, |a| {
            (PropValue::String(a.str(1)), 0.0)
        }),
        (
            "Double",
            "value +- tolerance",
            Kind::Double,
            Some(Param::optional(
                "tolerance",
                Kind::Double,
                self::Value::Double(0.0),
            )),
            |a| (PropValue::Double(a.double(1)), a.double(3)),
        ),
    ];
    for (name, doc, kind, tolerance, value) in kinds {
        for target in ["Atom", "Bond"] {
            let mut params = vec![
                Param::required("propname", Kind::Str),
                Param::required("val", kind),
                negate_param(),
            ];
            if let Some(tol) = &tolerance {
                params.push(Param {
                    name: tol.name,
                    kind: tol.kind,
                    default: tol.default.clone(),
                });
            }
            def_factory(
                m,
                format!("Has{name}PropWithValueQuery{target}"),
                format!(
                    "Returns a Query{target} that matches when the property 'propname' has the specified {doc}."
                ),
                params,
                move |py, a| {
                    let (val, tol) = value(a);
                    if target == "Atom" {
                        atom_object(
                            py,
                            make_prop_query::<AtomFn>(&a.str(0), val, tol),
                            a.bool(2),
                        )
                    } else {
                        bond_object(
                            py,
                            make_prop_query::<BondFn>(&a.str(0), val, tol),
                            a.bool(2),
                        )
                    }
                },
            )?;
        }
    }
    Ok(())
}

/// Changes the given atom in the molecule to
/// a query atom and returns the atom which can then be modified, for example
/// with additional query constraints added.  The new atom is otherwise a copy
/// of the old.
/// If the atom already has a query, nothing will be changed.
#[pyfunction]
#[pyo3(name = "ReplaceAtomWithQueryAtom", signature = (*args, **kwargs))]
fn py_replace_atom_with_query_atom(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Py<PyAtom>> {
    let call = Args {
        name: "ReplaceAtomWithQueryAtom",
        args,
        kwargs,
    };
    let mol = call.get(0, "mol")?;
    let atom = call.get(1, "atom")?;
    let (Some(Ok(mol)), Some(Ok(atom)), true) = (
        mol.as_ref().map(|m| m.cast::<PyMol>()),
        atom.as_ref().map(|a| a.cast::<PyAtom>()),
        args.len() + kwargs.map_or(0, |k| k.len()) == 2,
    ) else {
        let got: Vec<String> = args.iter().map(|a| boost_type_name(&a)).collect();
        return Err(argument_error(
            py,
            "ReplaceAtomWithQueryAtom",
            &got,
            "ReplaceAtomWithQueryAtom(RDKit::ROMol {lvalue} mol, RDKit::Atom {lvalue} atom)",
        ));
    };
    let (has_query, idx) = atom.borrow().with(py, |a| (a.has_query(), a.idx()))?;
    if has_query {
        return Ok(atom.clone().unbind());
    }
    replace_atom_with_query_atom(&mut mol.borrow_mut().mol, idx).map_err(crate::errors::to_py)?;
    PyAtom::in_mol(py, mol.clone().unbind(), idx)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing RDKit functionality for querying molecules.",
    )?;
    def_comparisons(m)?;
    def_flags(m)?;
    def_property_queries(m)?;
    m.add_function(wrap_pyfunction!(py_replace_atom_with_query_atom, m)?)?;
    Ok(())
}
