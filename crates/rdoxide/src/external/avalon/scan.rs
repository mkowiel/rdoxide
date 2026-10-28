//! Emulation of the `sscanf` conversions Avalon's MOL file reader relies on
//! (fixed-width `%Nd`, `%Nf`, `%Nc`, white space separated `%d`, `%f`, `%s`).
//! A conversion that fails stops the scan, as `sscanf` does.

/// A cursor over a byte string.
#[derive(Debug)]
pub struct Scan<'a> {
    s: &'a [u8],
    pub pos: usize,
}

impl<'a> Scan<'a> {
    pub fn new(s: &'a [u8]) -> Self {
        Scan { s, pos: 0 }
    }

    fn at(&self, i: usize) -> u8 {
        self.s.get(i).copied().unwrap_or(0)
    }

    fn skip_ws(&mut self) {
        while self.pos < self.s.len() && self.s[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }

    /// skips white space like a blank in the format
    pub fn blank(&mut self) {
        self.skip_ws();
    }

    /// `%d` (width 0: unlimited) after skipping white space
    pub fn int(&mut self, width: usize) -> Option<i32> {
        self.skip_ws();
        let limit = if width == 0 { usize::MAX } else { width };
        let mut i = self.pos;
        let mut used = 0;
        let neg = self.at(i) == b'-';
        if (self.at(i) == b'-' || self.at(i) == b'+') && used < limit {
            i += 1;
            used += 1;
        }
        let st = i;
        let mut v: i64 = 0;
        while self.at(i).is_ascii_digit() && used < limit {
            v = (v * 10 + i64::from(self.at(i) - b'0')).min(i64::from(u32::MAX));
            i += 1;
            used += 1;
        }
        if i == st {
            return None;
        }
        self.pos = i;
        Some(if neg { -(v as i32) } else { v as i32 })
    }

    /// `%f` (width 0: unlimited) after skipping white space
    pub fn float(&mut self, width: usize) -> Option<f32> {
        let (t, end) = self.float_token(width)?;
        let v = t.parse::<f32>().ok()?;
        self.pos = end;
        Some(v)
    }

    /// `%lf` (width 0: unlimited) after skipping white space
    pub fn double(&mut self, width: usize) -> Option<f64> {
        let (t, end) = self.float_token(width)?;
        let v = t.parse::<f64>().ok()?;
        self.pos = end;
        Some(v)
    }

    fn float_token(&mut self, width: usize) -> Option<(&'a str, usize)> {
        self.skip_ws();
        let limit = if width == 0 { usize::MAX } else { width };
        let start = self.pos;
        let mut i = start;
        let used = std::cell::Cell::new(0usize);
        let s = self.s;
        let at = |i: usize| s.get(i).copied().unwrap_or(0);
        let take = |i: &mut usize, pred: &dyn Fn(u8) -> bool| -> bool {
            if pred(at(*i)) && used.get() < limit {
                *i += 1;
                used.set(used.get() + 1);
                true
            } else {
                false
            }
        };
        take(&mut i, &|c| c == b'-' || c == b'+');
        let mut digits = 0;
        while take(&mut i, &|c| c.is_ascii_digit()) {
            digits += 1;
        }
        if take(&mut i, &|c| c == b'.') {
            while take(&mut i, &|c| c.is_ascii_digit()) {
                digits += 1;
            }
        }
        if digits == 0 {
            return None;
        }
        if at(i) == b'e' || at(i) == b'E' {
            let used_save = used.get();
            let mut j = i;
            if take(&mut j, &|_| true) {
                take(&mut j, &|c| c == b'-' || c == b'+');
                let mut ed = 0;
                while take(&mut j, &|c| c.is_ascii_digit()) {
                    ed += 1;
                }
                if ed > 0 {
                    i = j;
                } else {
                    used.set(used_save);
                }
            }
        }
        let t = std::str::from_utf8(&s[start..i]).ok()?;
        Some((t, i))
    }

    /// `%Nc`: exactly `n` characters
    pub fn chars(&mut self, n: usize) -> Option<String> {
        // a short read succeeds as long as one character is available
        if self.pos >= self.s.len() {
            return None;
        }
        let end = (self.pos + n).min(self.s.len());
        let t = String::from_utf8_lossy(&self.s[self.pos..end]).into_owned();
        self.pos = end;
        Some(t)
    }

    /// `%*Nc`: skips `n` characters
    pub fn skip(&mut self, n: usize) -> bool {
        if self.pos >= self.s.len() {
            return false;
        }
        self.pos = (self.pos + n).min(self.s.len());
        true
    }

    /// `%s`
    pub fn word(&mut self) -> Option<String> {
        self.skip_ws();
        let st = self.pos;
        while self.pos < self.s.len() && !self.s[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
        if st == self.pos {
            return None;
        }
        Some(String::from_utf8_lossy(&self.s[st..self.pos]).into_owned())
    }
}

/// C `printf("%g")` of a double
pub fn fmt_g(v: f64) -> String {
    fmt_gp(v, 6)
}

/// C `printf("%.<prec>g")` of a double
pub fn fmt_gp(v: f64, prec: usize) -> String {
    let prec = prec.max(1);
    if v == 0.0 {
        return if v.is_sign_negative() {
            "-0".into()
        } else {
            "0".into()
        };
    }
    if !v.is_finite() {
        return if v.is_nan() {
            "nan".into()
        } else if v < 0.0 {
            "-inf".into()
        } else {
            "inf".into()
        };
    }
    let e = format!("{:.*e}", prec - 1, v); // d.ddddde[-]x with prec significant digits
    let (mant, exp) = e.split_once('e').expect("exponent");
    let exp: i32 = exp.parse().expect("exponent");
    if exp < -4 || exp >= prec as i32 {
        let mant = if mant.contains('.') {
            mant.trim_end_matches('0').trim_end_matches('.')
        } else {
            mant
        };
        return format!(
            "{}e{}{:02}",
            mant,
            if exp < 0 { '-' } else { '+' },
            exp.abs()
        );
    }
    let dec = (prec as i32 - 1 - exp).max(0) as usize;
    let s = format!("{:.*}", dec, v);
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g_format() {
        assert_eq!(fmt_g(1.0), "1");
        assert_eq!(fmt_g(0.5), "0.5");
        assert_eq!(fmt_g(123456.0), "123456");
        assert_eq!(fmt_g(1234567.0), "1.23457e+06");
        assert_eq!(fmt_g(0.0001), "0.0001");
        assert_eq!(fmt_g(0.00001), "1e-05");
        assert_eq!(fmt_g(-2.25), "-2.25");
        assert_eq!(fmt_gp(12.5, 2), "12");
        assert_eq!(fmt_gp(99.7, 2), "1e+02");
        assert_eq!(fmt_gp(0.123, 2), "0.12");
    }

    #[test]
    fn fixed_width() {
        let mut s = Scan::new(b"  1  2 12");
        assert_eq!(s.int(3), Some(1));
        assert_eq!(s.int(3), Some(2));
        assert_eq!(s.int(3), Some(12));
        assert_eq!(s.int(3), None);
        let mut s = Scan::new(b"   -0.5000    1.2500");
        assert_eq!(s.float(10), Some(-0.5));
        assert_eq!(s.float(10), Some(1.25));
    }
}
