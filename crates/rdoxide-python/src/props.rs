//! Property access shared by Atom, Bond, Mol and Conformer — port of upstream
//! `Code/GraphMol/Wrap/props.hpp`.

use pyo3::BoundObject;
use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use rdoxide::rdgeneral::{Dict, RDProps, RDValue};

/// `from_rdvalue<T>` for arithmetic T: string values are lexically converted
/// (after trimming trailing whitespace).
fn lexical<T: std::str::FromStr>(v: &RDValue, direct: impl Fn(&RDValue) -> Option<T>) -> Option<T> {
    if let Some(x) = direct(v) {
        return Some(x);
    }
    if let RDValue::String(s) = v {
        return s.trim_end().parse::<T>().ok();
    }
    None
}

pub fn as_int(v: &RDValue) -> Option<i32> {
    lexical(v, RDValue::as_int)
}

pub fn as_uint(v: &RDValue) -> Option<u32> {
    // boost::lexical_cast<unsigned> accepts "-1" and wraps
    lexical(v, RDValue::as_uint).or_else(|| match v {
        RDValue::String(s) => s
            .trim_end()
            .parse::<i64>()
            .ok()
            .filter(|x| *x < 0 && *x >= -(u32::MAX as i64))
            .map(|x| x as u32),
        _ => None,
    })
}

/// `boost::lexical_cast<double>`: like Rust's parser, but an overflowing
/// exponent is a failure rather than infinity.
pub fn lexical_double(t: &str) -> Option<f64> {
    let d = t.parse::<f64>().ok()?;
    if d.is_infinite() {
        let u = t.trim_start_matches(['+', '-']).to_ascii_lowercase();
        if !u.starts_with("inf") {
            return None;
        }
    }
    Some(d)
}

pub fn as_double(v: &RDValue) -> Option<f64> {
    if let Some(x) = v.as_double() {
        return Some(x);
    }
    match v {
        RDValue::String(s) => lexical_double(s.trim_end()),
        _ => None,
    }
}

pub fn as_bool(v: &RDValue) -> Option<bool> {
    if let Some(b) = v.as_bool() {
        return Some(b);
    }
    match v {
        RDValue::String(s) => match s.trim_end() {
            "0" => Some(false),
            "1" => Some(true),
            _ => None,
        },
        _ => None,
    }
}

fn type_name(kind: &str) -> &'static str {
    match kind {
        "int" => "an integer value",
        "uint" => "an unsigned integer value",
        "double" => "a double value",
        "bool" => "a True or False value",
        _ => "unregistered C++ type",
    }
}

/// A property name; anything but `str` fails the way a Boost.Python overload mismatch does.
pub struct PropKey(pub String);

impl<'a, 'py> pyo3::FromPyObject<'a, 'py> for PropKey {
    type Error = PyErr;

    fn extract(obj: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        match obj.extract::<String>() {
            Ok(s) => Ok(Self(s)),
            Err(_) => Err(pyo3::exceptions::PyTypeError::new_err(format!(
                "Python argument types did not match C++ signature: property name must be str, got {}",
                crate::errors::boost_type_name(&obj)
            ))),
        }
    }
}

/// Why a stored number does not fit the requested integer type: boost's
/// `bad_numeric_cast` text for out-of-range numbers, `bad any cast` otherwise.
fn cast_failure_reason(v: &RDValue, kind: &str) -> &'static str {
    let wide = match v {
        RDValue::Int(i) => Some(i64::from(*i)),
        RDValue::UnsignedInt(u) => Some(i64::from(*u)),
        RDValue::Int64(i) => Some(*i),
        _ => None,
    };
    match (kind, wide) {
        ("uint", Some(x)) if x < 0 => "bad numeric conversion: negative overflow",
        ("int", Some(x)) if x > i64::from(i32::MAX) => "bad numeric conversion: positive overflow",
        _ => "bad any cast",
    }
}

/// Typed getter (`GetIntProp` & co.).
pub fn get_typed<'py, T: IntoPyObject<'py>>(
    py: Python<'py>,
    d: &Dict,
    key: &str,
    kind: &str,
    conv: fn(&RDValue) -> Option<T>,
    default: Option<Bound<'py, PyAny>>,
) -> PyResult<Py<PyAny>> {
    match d.get(key) {
        None => match default {
            Some(dv) => Ok(dv.unbind()),
            None => Err(PyKeyError::new_err(key.to_string())),
        },
        Some(v) => match conv(v) {
            Some(x) => Ok(x.into_pyobject(py).map_err(Into::into)?.into_any().unbind()),
            None => Err(PyValueError::new_err(format!(
                "key `{key}` exists but does not result in {} reason: {}",
                type_name(kind),
                cast_failure_reason(v, kind)
            ))),
        },
    }
}

/// Python object for a stored value (autoConvert / GetPropsAsDict).
pub fn rdvalue_to_py(
    py: Python<'_>,
    v: &RDValue,
    auto_convert_strings: bool,
) -> PyResult<Py<PyAny>> {
    Ok(match v {
        RDValue::BitVect(bv) => {
            Py::new(py, crate::datastructs::PyExplicitBitVect::wrap(bv.clone()))?.into_any()
        }
        RDValue::Int(i) => i.into_pyobject(py)?.into_any().unbind(),
        RDValue::UnsignedInt(u) => u.into_pyobject(py)?.into_any().unbind(),
        RDValue::Int64(i) => i.into_pyobject(py)?.into_any().unbind(),
        RDValue::UInt64(u) => u.into_pyobject(py)?.into_any().unbind(),
        RDValue::Double(d) => d.into_pyobject(py)?.into_any().unbind(),
        RDValue::Float(f) => f.into_pyobject(py)?.into_any().unbind(),
        RDValue::Bool(b) => b.into_pyobject(py)?.to_owned().into_any().unbind(),
        RDValue::String(s) => {
            if auto_convert_strings {
                let t = s.trim();
                if let Ok(i) = t.parse::<i32>() {
                    return Ok(i.into_pyobject(py)?.into_any().unbind());
                }
                if let Some(d) = lexical_double(t) {
                    return Ok(d.into_pyobject(py)?.into_any().unbind());
                }
            }
            s.into_pyobject(py)?.into_any().unbind()
        }
        RDValue::VecInt(x) => crate::vects::PyVectI::new_py(py, x.clone())?,
        RDValue::VecUnsignedInt(x) => x.clone().into_pyobject(py)?.into_any().unbind(),
        RDValue::VecDouble(x) => x.clone().into_pyobject(py)?.into_any().unbind(),
        RDValue::VecFloat(x) => x.clone().into_pyobject(py)?.into_any().unbind(),
        RDValue::VecString(x) => crate::vects::PyVectStr::new_py(py, x.clone())?,
    })
}

/// `autoConvertString`: int, then double, then the string itself (no trim).
fn auto_convert_string(py: Python<'_>, v: &RDValue) -> PyResult<Py<PyAny>> {
    if let Some(i) = as_int(v) {
        return Ok(i.into_pyobject(py)?.into_any().unbind());
    }
    if let Some(d) = as_double(v) {
        return Ok(d.into_pyobject(py)?.into_any().unbind());
    }
    Ok(v.to_rd_string().into_pyobject(py)?.into_any().unbind())
}

/// `GetPyPropImpl`
pub fn get_py_prop(
    py: Python<'_>,
    d: &Dict,
    key: &str,
    auto_convert: bool,
    default: Option<Bound<'_, PyAny>>,
) -> PyResult<Py<PyAny>> {
    match d.get(key) {
        None => match default {
            Some(dv) => Ok(dv.unbind()),
            None => Err(PyKeyError::new_err(key.to_string())),
        },
        Some(v) if !auto_convert => Ok(v.to_rd_string().into_pyobject(py)?.into_any().unbind()),
        Some(v @ RDValue::String(_)) => auto_convert_string(py, v),
        Some(v) => rdvalue_to_py(py, v, false),
    }
}

/// `GetPropsAsDict`
pub fn props_as_dict<'py>(
    py: Python<'py>,
    d: &Dict,
    include_private: bool,
    include_computed: bool,
    auto_convert_strings: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let res = PyDict::new(py);
    let keys = d.get_prop_list(include_private, include_computed);
    for (k, v) in d.iter() {
        if keys.iter().any(|x| x == k) && !rdoxide::graphmol::matrices::is_matrix_cache_prop(k) {
            res.set_item(k, rdvalue_to_py(py, v, auto_convert_strings)?)?;
        }
    }
    Ok(res)
}

/// Fetch an optional argument by position (after the fixed ones) or by
/// keyword, distinguishing "not given" from an explicit `None`.
pub fn opt_arg<'py>(
    args: &Bound<'py, pyo3::types::PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
    pos: usize,
    name: &str,
) -> PyResult<Option<Bound<'py, PyAny>>> {
    if args.len() > pos {
        return Ok(Some(args.get_item(pos)?));
    }
    match kwargs {
        Some(kw) => kw.get_item(name),
        None => Ok(None),
    }
}

pub fn set_prop(d: &mut Dict, key: &str, val: RDValue, computed: bool) -> PyResult<()> {
    d.set_prop(key, val, computed).map_err(crate::errors::to_py)
}

/// Generates the property methods for a pyclass that implements
/// `with_dict(&self, py, |&Dict| ..)` and `with_dict_mut(&self, py, |&mut Dict| ..)`.
///
/// `computed`: whether setters take a `computed` argument (Mol, Conformer).
/// `$priv`, `$comp`: defaults of `GetPropsAsDict(includePrivate, includeComputed)`.
/// `has_prop`: the type `HasProp` returns (`i32` for most classes, `bool` where upstream does).
macro_rules! impl_props {
    ($ty:ty, computed = $computed:tt, dict_defaults = ($priv:expr, $comp:expr)) => {
        $crate::props::impl_props!($ty, computed = $computed, dict_defaults = ($priv, $comp), has_prop = i32);
    };
    ($ty:ty, computed = $computed:tt, dict_defaults = ($priv:expr, $comp:expr), has_prop = $has:ty) => {
        $crate::props::impl_props!(@setters $ty, $computed);

        #[pyo3::pymethods]
        impl $ty {
            /// Queries a molecule to see if a particular property has been assigned.
            #[pyo3(name = "HasProp")]
            fn py_has_prop(&self, py: pyo3::Python<'_>, key: $crate::props::PropKey) -> pyo3::PyResult<$has> {
                self.with_dict(py, |d| <$has>::from(d.has_val(&key.0)))
            }

            /// Removes a particular property (does nothing if not already set).
            #[pyo3(name = "ClearProp")]
            fn py_clear_prop(&mut self, py: pyo3::Python<'_>, key: $crate::props::PropKey) -> pyo3::PyResult<()> {
                self.with_dict_mut(py, |d| {
                    if d.has_val(&key.0) {
                        rdoxide::rdgeneral::RDProps::clear_prop(d, &key.0);
                    }
                })
            }

            /// Returns the value of the property.
            #[pyo3(name = "GetProp", signature = (key, *args, **kwargs))]
            fn py_get_prop(
                &self,
                py: pyo3::Python<'_>,
                key: $crate::props::PropKey,
                args: &pyo3::Bound<'_, pyo3::types::PyTuple>,
                kwargs: Option<&pyo3::Bound<'_, pyo3::types::PyDict>>,
            ) -> pyo3::PyResult<pyo3::Py<pyo3::PyAny>> {
                use pyo3::types::PyAnyMethods;
                let auto_convert = match $crate::props::opt_arg(args, kwargs, 0, "autoConvert")? {
                    Some(v) => v.is_truthy()?,
                    None => false,
                };
                let default = $crate::props::opt_arg(args, kwargs, 1, "default")?;
                self.with_dict(py, |d| $crate::props::get_py_prop(py, d, &key.0, auto_convert, default))?
            }

            /// Returns the integer value of the property if possible.
            #[pyo3(name = "GetIntProp", signature = (key, *args, **kwargs))]
            fn py_get_int_prop(&self, py: pyo3::Python<'_>, key: $crate::props::PropKey, args: &pyo3::Bound<'_, pyo3::types::PyTuple>, kwargs: Option<&pyo3::Bound<'_, pyo3::types::PyDict>>) -> pyo3::PyResult<pyo3::Py<pyo3::PyAny>> {
                let default = $crate::props::opt_arg(args, kwargs, 0, "default")?;
                self.with_dict(py, |d| $crate::props::get_typed(py, d, &key.0, "int", $crate::props::as_int, default))?
            }

            /// Returns the unsigned int value of the property if possible.
            #[pyo3(name = "GetUnsignedProp", signature = (key, *args, **kwargs))]
            fn py_get_unsigned_prop(&self, py: pyo3::Python<'_>, key: $crate::props::PropKey, args: &pyo3::Bound<'_, pyo3::types::PyTuple>, kwargs: Option<&pyo3::Bound<'_, pyo3::types::PyDict>>) -> pyo3::PyResult<pyo3::Py<pyo3::PyAny>> {
                let default = $crate::props::opt_arg(args, kwargs, 0, "default")?;
                self.with_dict(py, |d| $crate::props::get_typed(py, d, &key.0, "uint", $crate::props::as_uint, default))?
            }

            /// Returns the double value of the property if possible.
            #[pyo3(name = "GetDoubleProp", signature = (key, *args, **kwargs))]
            fn py_get_double_prop(&self, py: pyo3::Python<'_>, key: $crate::props::PropKey, args: &pyo3::Bound<'_, pyo3::types::PyTuple>, kwargs: Option<&pyo3::Bound<'_, pyo3::types::PyDict>>) -> pyo3::PyResult<pyo3::Py<pyo3::PyAny>> {
                let default = $crate::props::opt_arg(args, kwargs, 0, "default")?;
                self.with_dict(py, |d| $crate::props::get_typed(py, d, &key.0, "double", $crate::props::as_double, default))?
            }

            /// Returns the Bool value of the property if possible.
            #[pyo3(name = "GetBoolProp", signature = (key, *args, **kwargs))]
            fn py_get_bool_prop(&self, py: pyo3::Python<'_>, key: $crate::props::PropKey, args: &pyo3::Bound<'_, pyo3::types::PyTuple>, kwargs: Option<&pyo3::Bound<'_, pyo3::types::PyDict>>) -> pyo3::PyResult<pyo3::Py<pyo3::PyAny>> {
                let default = $crate::props::opt_arg(args, kwargs, 0, "default")?;
                self.with_dict(py, |d| $crate::props::get_typed(py, d, &key.0, "bool", $crate::props::as_bool, default))?
            }

            /// Returns a list of the properties set.
            #[pyo3(name = "GetPropNames", signature = (includePrivate = false, includeComputed = false))]
            #[allow(non_snake_case)]
            fn py_get_prop_names(&self, py: pyo3::Python<'_>, includePrivate: bool, includeComputed: bool) -> pyo3::PyResult<$crate::vects::PyVectStr> {
                let v = self.with_dict(py, |d| rdoxide::rdgeneral::RDProps::get_prop_list(d, includePrivate, includeComputed))?;
                Ok($crate::vects::PyVectStr { v })
            }

            /// Returns a dictionary populated with properties.
            /// When possible, string values will be converted to integers or doubles (trimming if necessary)
            #[pyo3(name = "GetPropsAsDict", signature = (includePrivate = $priv, includeComputed = $comp, autoConvertStrings = true))]
            #[allow(non_snake_case)]
            fn py_get_props_as_dict<'py>(
                &self,
                py: pyo3::Python<'py>,
                includePrivate: bool,
                includeComputed: bool,
                autoConvertStrings: bool,
            ) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyDict>> {
                self.with_dict(py, |d| $crate::props::props_as_dict(py, d, includePrivate, includeComputed, autoConvertStrings))?
            }
        }
    };
    (@setters $ty:ty, true) => {
        #[pyo3::pymethods]
        impl $ty {
            /// Sets a property (string value).
            #[pyo3(name = "SetProp", signature = (key, val, computed = false))]
            fn py_set_prop(&mut self, py: pyo3::Python<'_>, key: $crate::props::PropKey, val: String, #[pyo3(from_py_with = crate::pyconv::truthy)] computed: bool) -> pyo3::PyResult<()> {
                self.with_dict_mut(py, |d| $crate::props::set_prop(d, &key.0, val.into(), computed))?
            }
            /// Sets a double valued property.
            #[pyo3(name = "SetDoubleProp", signature = (key, val, computed = false))]
            fn py_set_double_prop(&mut self, py: pyo3::Python<'_>, key: $crate::props::PropKey, val: f64, #[pyo3(from_py_with = crate::pyconv::truthy)] computed: bool) -> pyo3::PyResult<()> {
                self.with_dict_mut(py, |d| $crate::props::set_prop(d, &key.0, val.into(), computed))?
            }
            /// Sets an integer valued property.
            #[pyo3(name = "SetIntProp", signature = (key, val, computed = false))]
            fn py_set_int_prop(&mut self, py: pyo3::Python<'_>, key: $crate::props::PropKey, val: i32, #[pyo3(from_py_with = crate::pyconv::truthy)] computed: bool) -> pyo3::PyResult<()> {
                self.with_dict_mut(py, |d| $crate::props::set_prop(d, &key.0, val.into(), computed))?
            }
            /// Sets an unsigned integer valued property.
            #[pyo3(name = "SetUnsignedProp", signature = (key, val, computed = false))]
            fn py_set_unsigned_prop(&mut self, py: pyo3::Python<'_>, key: $crate::props::PropKey, val: u32, #[pyo3(from_py_with = crate::pyconv::truthy)] computed: bool) -> pyo3::PyResult<()> {
                self.with_dict_mut(py, |d| $crate::props::set_prop(d, &key.0, val.into(), computed))?
            }
            /// Sets a boolean valued property.
            #[pyo3(name = "SetBoolProp", signature = (key, val, computed = false))]
            fn py_set_bool_prop(&mut self, py: pyo3::Python<'_>, key: $crate::props::PropKey, #[pyo3(from_py_with = crate::pyconv::truthy)] val: bool, #[pyo3(from_py_with = crate::pyconv::truthy)] computed: bool) -> pyo3::PyResult<()> {
                self.with_dict_mut(py, |d| $crate::props::set_prop(d, &key.0, val.into(), computed))?
            }
        }
    };
    (@setters $ty:ty, false) => {
        #[pyo3::pymethods]
        impl $ty {
            /// Sets a property (string value).
            #[pyo3(name = "SetProp")]
            fn py_set_prop(&mut self, py: pyo3::Python<'_>, key: $crate::props::PropKey, val: String) -> pyo3::PyResult<()> {
                self.with_dict_mut(py, |d| $crate::props::set_prop(d, &key.0, val.into(), false))?
            }
            /// Sets a double valued property.
            #[pyo3(name = "SetDoubleProp")]
            fn py_set_double_prop(&mut self, py: pyo3::Python<'_>, key: $crate::props::PropKey, val: f64) -> pyo3::PyResult<()> {
                self.with_dict_mut(py, |d| $crate::props::set_prop(d, &key.0, val.into(), false))?
            }
            /// Sets an integer valued property.
            #[pyo3(name = "SetIntProp")]
            fn py_set_int_prop(&mut self, py: pyo3::Python<'_>, key: $crate::props::PropKey, val: i32) -> pyo3::PyResult<()> {
                self.with_dict_mut(py, |d| $crate::props::set_prop(d, &key.0, val.into(), false))?
            }
            /// Sets an unsigned integer valued property.
            #[pyo3(name = "SetUnsignedProp")]
            fn py_set_unsigned_prop(&mut self, py: pyo3::Python<'_>, key: $crate::props::PropKey, val: u32) -> pyo3::PyResult<()> {
                self.with_dict_mut(py, |d| $crate::props::set_prop(d, &key.0, val.into(), false))?
            }
            /// Sets a boolean valued property.
            #[pyo3(name = "SetBoolProp")]
            fn py_set_bool_prop(&mut self, py: pyo3::Python<'_>, key: $crate::props::PropKey, val: bool) -> pyo3::PyResult<()> {
                self.with_dict_mut(py, |d| $crate::props::set_prop(d, &key.0, val.into(), false))?
            }
        }
    };
}
pub(crate) use impl_props;
