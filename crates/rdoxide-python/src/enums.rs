//! Boost.Python-compatible enum types.
//!
//! RDKit's released wheels expose C++ enums through Boost.Python: the enum
//! class is an `int` subclass, `str(x)` is the value's name, `repr(x)` is
//! `module.Enum.NAME`, and the class has `names`/`values` dictionaries. That
//! behaviour is reproduced by a tiny Python factory evaluated at import time;
//! Rust code converts values through the per-enum registry below.

use std::collections::HashMap;
use std::ffi::CString;
use std::sync::Mutex;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};

const FACTORY: &str = r#"
def make_enum(module, name, members, doc=None):
    def __new__(cls, value=0):
        v = int(value)
        try:
            return cls.values[v]
        except KeyError:
            obj = int.__new__(cls, v)
            obj.name = None
            return obj

    def __repr__(self):
        if self.name is None:
            return "%s.%s(%d)" % (module, name, int(self))
        return "%s.%s.%s" % (module, name, self.name)

    def __str__(self):
        return self.name if self.name is not None else str(int(self))

    def __reduce__(self):
        return (type(self), (int(self),))

    ns = {
        "__new__": __new__,
        "__repr__": __repr__,
        "__str__": __str__,
        "__reduce__": __reduce__,
        "__module__": module,
        "__doc__": doc,
        "__hash__": int.__hash__,
        "__eq__": int.__eq__,
        "__ne__": int.__ne__,
    }
    cls = type(name, (int,), ns)
    cls.names = {}
    cls.values = {}
    for n, v in members:
        obj = int.__new__(cls, v)
        obj.name = n
        setattr(cls, n, obj)
        cls.names[n] = obj
        cls.values.setdefault(v, obj)
    return cls
"#;

struct EnumEntry {
    cls: Py<PyAny>,
    values: HashMap<i64, Py<PyAny>>,
}

static REGISTRY: Mutex<Option<HashMap<&'static str, EnumEntry>>> = Mutex::new(None);

/// Create the enum class `name` in module `m` and register it.
pub fn add_enum(
    m: &Bound<'_, PyModule>,
    name: &'static str,
    members: &[(&str, i64)],
    doc: Option<&str>,
) -> PyResult<()> {
    let py = m.py();
    let code = CString::new(FACTORY).expect("no NUL in factory code");
    let helper = PyModule::from_code(py, code.as_c_str(), c"rdoxide_enum.py", c"rdoxide_enum")?;
    let module_name: String = m.getattr("__name__")?.extract()?;
    let items = PyList::empty(py);
    for (n, v) in members {
        items.append(PyTuple::new(
            py,
            [
                n.into_pyobject(py)?.into_any(),
                v.into_pyobject(py)?.into_any(),
            ],
        )?)?;
    }
    let cls = helper
        .getattr("make_enum")?
        .call1((module_name, name, items, doc))?;
    let values_dict = cls.getattr("values")?.cast_into::<PyDict>()?;
    let mut values = HashMap::new();
    for (k, v) in values_dict.iter() {
        values.insert(k.extract::<i64>()?, v.unbind());
    }
    m.add(name, &cls)?;
    let mut reg = REGISTRY.lock().unwrap_or_else(|e| e.into_inner());
    reg.get_or_insert_with(HashMap::new).insert(
        name,
        EnumEntry {
            cls: cls.unbind(),
            values,
        },
    );
    Ok(())
}

/// Like [`add_enum`] and also export the members into the module namespace
/// (Boost.Python's `.export_values()`).
pub fn add_enum_exported(
    m: &Bound<'_, PyModule>,
    name: &'static str,
    members: &[(&str, i64)],
    doc: Option<&str>,
) -> PyResult<()> {
    add_enum(m, name, members, doc)?;
    let cls = m.getattr(name)?;
    for (n, _) in members {
        m.add(*n, cls.getattr(*n)?)?;
    }
    Ok(())
}

/// Convert a value to an instance of the registered enum `name`.
pub fn to_py(py: Python<'_>, name: &str, value: i64) -> PyResult<Py<PyAny>> {
    let cls = {
        let reg = REGISTRY.lock().unwrap_or_else(|e| e.into_inner());
        let entry = reg
            .as_ref()
            .and_then(|r| r.get(name))
            .ok_or_else(|| PyValueError::new_err(format!("unknown enum {name}")))?;
        if let Some(v) = entry.values.get(&value) {
            return Ok(v.clone_ref(py));
        }
        entry.cls.clone_ref(py)
    };
    Ok(cls.bind(py).call1((value,))?.unbind())
}

/// Declare a Rust enum <-> Python enum mapping.
macro_rules! py_enum {
    ($wrapper:ident, $rust:ty, $pyname:literal) => {
        #[derive(Debug, Clone, Copy)]
        pub struct $wrapper(pub $rust);

        impl<'py> pyo3::IntoPyObject<'py> for $wrapper {
            type Target = pyo3::PyAny;
            type Output = pyo3::Bound<'py, pyo3::PyAny>;
            type Error = pyo3::PyErr;
            fn into_pyobject(self, py: pyo3::Python<'py>) -> Result<Self::Output, Self::Error> {
                Ok($crate::enums::to_py(py, $pyname, self.0 as i64)?.into_bound(py))
            }
        }

        impl $wrapper {
            /// Extract from any Python int (including enum instances).
            pub fn extract(obj: &pyo3::Bound<'_, pyo3::PyAny>) -> pyo3::PyResult<$rust> {
                let v: i64 = pyo3::types::PyAnyMethods::extract(obj)?;
                u8::try_from(v)
                    .ok()
                    .and_then(<$rust>::from_u8)
                    .ok_or_else(|| {
                        pyo3::exceptions::PyValueError::new_err(format!(
                            "invalid {} value: {}",
                            $pyname, v
                        ))
                    })
            }
            pub fn members() -> Vec<(&'static str, i64)> {
                <$rust>::ALL.iter().map(|v| (v.name(), *v as i64)).collect()
            }
        }
    };
}
pub(crate) use py_enum;
