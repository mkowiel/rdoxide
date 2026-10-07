//! Ported from RDKit `Code/RDGeneral/Dict.h`, `RDValue*.h` and `RDProps.h`.
//!
//! Copyright (C) 2003-2024 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::fmt;

use super::strings::{double_to_string, float_to_string, format_g};

/// Name of the property holding the list of computed property names.
pub const COMPUTED_PROP_NAME: &str = "__computedProps";

/// A dynamically typed property value (`RDValue`).
#[derive(Debug, Clone, PartialEq)]
pub enum RDValue {
    Int(i32),
    UnsignedInt(u32),
    Int64(i64),
    UInt64(u64),
    Double(f64),
    Float(f32),
    Bool(bool),
    String(String),
    VecInt(Vec<i32>),
    VecUnsignedInt(Vec<u32>),
    VecDouble(Vec<f64>),
    VecFloat(Vec<f32>),
    VecString(Vec<String>),
}

/// Error returned when a value cannot be converted to the requested type
/// (`boost::bad_any_cast` / `bad_lexical_cast` upstream).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BadCast {
    pub key: String,
    pub requested: &'static str,
}

impl fmt::Display for BadCast {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "key `{}` exists but does not result in {} value",
            self.key, self.requested
        )
    }
}

impl std::error::Error for BadCast {}

fn vect_to_string<T>(v: &[T], fmt_one: impl Fn(&T) -> String) -> String {
    let parts: Vec<String> = v.iter().map(fmt_one).collect();
    format!("[{}]", parts.join(","))
}

impl RDValue {
    /// `rdvalue_tostring`: string representation used by `GetProp`.
    pub fn to_rd_string(&self) -> String {
        match self {
            RDValue::String(s) => s.clone(),
            RDValue::Int(i) => i.to_string(),
            RDValue::UnsignedInt(u) => u.to_string(),
            RDValue::Int64(i) => i.to_string(),
            RDValue::UInt64(u) => u.to_string(),
            RDValue::Double(d) => double_to_string(*d),
            RDValue::Float(f) => float_to_string(*f),
            RDValue::Bool(b) => if *b { "1" } else { "0" }.to_string(),
            RDValue::VecDouble(v) => vect_to_string(v, |x| format_g(*x, 17)),
            RDValue::VecFloat(v) => vect_to_string(v, |x| format_g(f64::from(*x), 17)),
            RDValue::VecInt(v) => vect_to_string(v, |x| x.to_string()),
            RDValue::VecUnsignedInt(v) => vect_to_string(v, |x| x.to_string()),
            RDValue::VecString(v) => vect_to_string(v, |x| x.clone()),
        }
    }

    pub fn as_int(&self) -> Option<i32> {
        match *self {
            RDValue::Int(i) => Some(i),
            RDValue::UnsignedInt(u) => i32::try_from(u).ok(),
            RDValue::Int64(i) => i32::try_from(i).ok(),
            RDValue::UInt64(u) => i32::try_from(u).ok(),
            _ => None,
        }
    }

    pub fn as_uint(&self) -> Option<u32> {
        match *self {
            RDValue::UnsignedInt(u) => Some(u),
            RDValue::Int(i) => u32::try_from(i).ok(),
            RDValue::Int64(i) => u32::try_from(i).ok(),
            RDValue::UInt64(u) => u32::try_from(u).ok(),
            _ => None,
        }
    }

    pub fn as_double(&self) -> Option<f64> {
        match *self {
            RDValue::Double(d) => Some(d),
            RDValue::Float(f) => Some(f64::from(f)),
            _ => None,
        }
    }

    pub fn as_bool(&self) -> Option<bool> {
        match *self {
            RDValue::Bool(b) => Some(b),
            _ => None,
        }
    }
}

macro_rules! impl_from {
    ($t:ty, $v:ident) => {
        impl From<$t> for RDValue {
            fn from(x: $t) -> Self {
                RDValue::$v(x)
            }
        }
    };
}
impl_from!(i32, Int);
impl_from!(u32, UnsignedInt);
impl_from!(i64, Int64);
impl_from!(u64, UInt64);
impl_from!(f64, Double);
impl_from!(f32, Float);
impl_from!(bool, Bool);
impl_from!(String, String);
impl_from!(Vec<i32>, VecInt);
impl_from!(Vec<u32>, VecUnsignedInt);
impl_from!(Vec<f64>, VecDouble);
impl_from!(Vec<f32>, VecFloat);
impl_from!(Vec<String>, VecString);

impl From<&str> for RDValue {
    fn from(x: &str) -> Self {
        RDValue::String(x.to_string())
    }
}

/// Insertion-ordered key/value store (`RDKit::Dict`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Dict {
    data: Vec<(String, RDValue)>,
}

impl Dict {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn keys(&self) -> impl Iterator<Item = &str> {
        self.data.iter().map(|(k, _)| k.as_str())
    }

    pub fn iter(&self) -> impl Iterator<Item = (&str, &RDValue)> {
        self.data.iter().map(|(k, v)| (k.as_str(), v))
    }

    pub fn has_val(&self, key: &str) -> bool {
        self.data.iter().any(|(k, _)| k == key)
    }

    pub fn get(&self, key: &str) -> Option<&RDValue> {
        self.data.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    pub fn get_mut(&mut self, key: &str) -> Option<&mut RDValue> {
        self.data.iter_mut().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    pub fn set_val(&mut self, key: &str, val: impl Into<RDValue>) {
        let val = val.into();
        match self.get_mut(key) {
            Some(slot) => *slot = val,
            None => self.data.push((key.to_string(), val)),
        }
    }

    pub fn clear_val(&mut self, key: &str) {
        self.data.retain(|(k, _)| k != key);
    }

    pub fn reset(&mut self) {
        self.data.clear();
    }

    /// `Dict::update`
    pub fn update(&mut self, other: &Dict, preserve_existing: bool) {
        for (k, v) in &other.data {
            if preserve_existing && self.has_val(k) {
                continue;
            }
            self.set_val(k, v.clone());
        }
    }
}

/// Property semantics shared by atoms, bonds, molecules, conformers
/// (`RDKit::RDProps`).
pub trait RDProps {
    fn dict(&self) -> &Dict;
    fn dict_mut(&mut self) -> &mut Dict;

    fn get_prop_list(&self, include_private: bool, include_computed: bool) -> Vec<String> {
        let mut computed: Vec<String> = Vec::new();
        if !include_computed {
            if let Some(RDValue::VecString(c)) = self.dict().get(COMPUTED_PROP_NAME) {
                computed = c.clone();
            }
            computed.push(COMPUTED_PROP_NAME.to_string());
        }
        self.dict()
            .keys()
            .filter(|k| {
                (include_private || !k.starts_with('_')) && !computed.iter().any(|c| c == k)
            })
            .map(str::to_string)
            .collect()
    }

    fn set_prop(
        &mut self,
        key: &str,
        val: impl Into<RDValue>,
        computed: bool,
    ) -> Result<(), crate::RDKitError>
    where
        Self: Sized,
    {
        if key.is_empty() {
            return Err(crate::RDKitError::Value(
                "Cannot set property with empty key".into(),
            ));
        }
        if computed {
            let mut lst = match self.dict().get(COMPUTED_PROP_NAME) {
                Some(RDValue::VecString(c)) => c.clone(),
                _ => Vec::new(),
            };
            if !lst.iter().any(|c| c == key) {
                lst.push(key.to_string());
                self.dict_mut().set_val(COMPUTED_PROP_NAME, lst);
            }
        }
        self.dict_mut().set_val(key, val);
        Ok(())
    }

    fn has_prop(&self, key: &str) -> bool {
        self.dict().has_val(key)
    }

    fn get_prop(&self, key: &str) -> Option<&RDValue> {
        self.dict().get(key)
    }

    fn clear_prop(&mut self, key: &str) {
        if let Some(RDValue::VecString(c)) = self.dict_mut().get_mut(COMPUTED_PROP_NAME) {
            c.retain(|x| x != key);
        }
        self.dict_mut().clear_val(key);
    }

    fn clear_computed_props(&mut self) {
        let lst = match self.dict().get(COMPUTED_PROP_NAME) {
            Some(RDValue::VecString(c)) if !c.is_empty() => c.clone(),
            _ => return,
        };
        for k in &lst {
            self.dict_mut().clear_val(k);
        }
        self.dict_mut()
            .set_val(COMPUTED_PROP_NAME, Vec::<String>::new());
    }

    fn update_props(&mut self, source: &Dict, preserve_existing: bool) {
        self.dict_mut().update(source, preserve_existing);
    }
}

impl RDProps for Dict {
    fn dict(&self) -> &Dict {
        self
    }
    fn dict_mut(&mut self) -> &mut Dict {
        self
    }
}

#[cfg(test)]
mod tests {
    //! Ported from Code/RDGeneral/testDict.cpp and catch_dict.cpp (subset).
    use super::*;

    #[test]
    fn basic_set_get() {
        let mut d = Dict::new();
        d.set_val("foo", 1);
        d.set_val("bar", 3.5);
        d.set_val("baz", "text");
        assert!(d.has_val("foo"));
        assert_eq!(d.get("foo").and_then(RDValue::as_int), Some(1));
        assert_eq!(d.get("bar").and_then(RDValue::as_double), Some(3.5));
        assert_eq!(d.get("baz").unwrap().to_rd_string(), "text");
        d.set_val("foo", 2);
        assert_eq!(d.len(), 3);
        assert_eq!(d.keys().collect::<Vec<_>>(), ["foo", "bar", "baz"]);
        d.clear_val("bar");
        assert!(!d.has_val("bar"));
    }

    #[test]
    fn string_conversions() {
        assert_eq!(RDValue::Int(-3).to_rd_string(), "-3");
        assert_eq!(RDValue::Double(0.1).to_rd_string(), "0.10000000000000001");
        assert_eq!(RDValue::Bool(true).to_rd_string(), "1");
        assert_eq!(RDValue::VecInt(vec![1, 2, 3]).to_rd_string(), "[1,2,3]");
        assert_eq!(RDValue::VecDouble(vec![]).to_rd_string(), "[]");
        assert_eq!(
            RDValue::VecString(vec!["a".into(), "b".into()]).to_rd_string(),
            "[a,b]"
        );
    }

    #[test]
    fn computed_props() {
        let mut d = Dict::new();
        d.set_prop("a", 1, false).unwrap();
        d.set_prop("b", 2, true).unwrap();
        d.set_prop("_c", 3, false).unwrap();
        assert_eq!(d.get_prop_list(false, false), ["a"]);
        assert_eq!(d.get_prop_list(true, false), ["a", "_c"]);
        assert_eq!(d.get_prop_list(false, true), ["a", "b"]);
        d.clear_computed_props();
        assert!(!d.has_prop("b"));
        assert!(d.set_prop("", 1, false).is_err());
    }
}
