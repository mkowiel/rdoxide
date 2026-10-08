//! Floating-point contraction as the reference C++ builds do it.
//!
//! Clang's default `-ffp-contract=on` fuses `a * b + c` into one
//! fused multiply-add when both appear in the same C++ expression, and on
//! targets with FMA instructions (arm64) this changes the last bit of
//! results. RDKit's arm64 builds (e.g. conda-forge macOS) therefore compute
//! `x * y + z` as `fma(x, y, z)`, while x86-64 builds (no FMA in the
//! baseline instruction set) round twice. [`fmuladd`] reproduces the
//! native build of the target platform; call sites follow clang's choice
//! of operands: in `a * b + c * d` the left product is fused, a lone
//! right-hand product is fused when the left operand is not a product, and
//! `a * b - c * d` becomes `fma(a, b, -(c * d))`.

/// `llvm.fmuladd(a, b, c)`: `a * b + c`, fused on aarch64.
#[inline(always)]
pub fn fmuladd(a: f64, b: f64, c: f64) -> f64 {
    #[cfg(target_arch = "aarch64")]
    {
        a.mul_add(b, c)
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        a * b + c
    }
}
