//! Number formatting compatible with C++ iostreams / `boost::lexical_cast`.
//!
//! RDKit converts numeric properties to strings with `boost::lexical_cast`
//! (which uses `max_digits10` significant digits, i.e. `%.17g` for doubles and
//! `%.9g` for floats) and with `std::ostream << std::setprecision(n)`
//! (`%.{n}g`). Rust has no `%g`, so it is implemented here.

/// Format `x` like C's `printf("%.{precision}g", x)`.
pub fn format_g(x: f64, precision: usize) -> String {
    if x.is_nan() {
        return if x.is_sign_negative() { "-nan" } else { "nan" }.to_string();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf" } else { "inf" }.to_string();
    }
    let p = precision.max(1);
    if x == 0.0 {
        return if x.is_sign_negative() { "-0" } else { "0" }.to_string();
    }
    // Determine the decimal exponent after rounding to p significant digits.
    let sci = format!("{:.*e}", p - 1, x);
    let (mantissa, exp) = sci.split_once('e').expect("exponent present");
    let exp: i32 = exp.parse().expect("valid exponent");
    if exp < -4 || exp >= p as i32 {
        let mantissa = strip_trailing_zeros(mantissa);
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{mantissa}e{sign}{:02}", exp.abs())
    } else {
        let decimals = (p as i32 - 1 - exp).max(0) as usize;
        strip_trailing_zeros(&format!("{x:.decimals$}")).to_string()
    }
}

fn strip_trailing_zeros(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

/// `boost::lexical_cast<std::string>(double)`
pub fn double_to_string(x: f64) -> String {
    format_g(x, 17)
}

/// `boost::lexical_cast<std::string>(float)`
pub fn float_to_string(x: f32) -> String {
    format_g(f64::from(x), 9)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_printf_g() {
        assert_eq!(format_g(1.5, 17), "1.5");
        assert_eq!(format_g(0.1, 17), "0.10000000000000001");
        assert_eq!(format_g(100.0, 17), "100");
        assert_eq!(format_g(1e20, 17), "1e+20");
        assert_eq!(format_g(1.234e-5, 17), "1.234e-05");
        assert_eq!(format_g(0.0001, 6), "0.0001");
        assert_eq!(format_g(123456.0, 6), "123456");
        assert_eq!(format_g(1234567.0, 6), "1.23457e+06");
        assert_eq!(format_g(-2.5, 6), "-2.5");
        assert_eq!(format_g(0.0, 17), "0");
        assert_eq!(float_to_string(0.1), "0.100000001");
        assert_eq!(format_g(9.9999999, 3), "10");
    }
}
