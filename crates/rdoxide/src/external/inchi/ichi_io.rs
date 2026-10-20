//! InChI I/O streams (`INCHI_BASE/src/ichi_io.{h,c}`): the resizable string
//! buffer (`INCHI_IOS_STRING`), the stream that is either such a buffer or a
//! plain file (`INCHI_IOSTREAM`), line readers and the growing print buffer.
//!
//! Differences from the C code: `printf`-style formatting is done by the
//! caller (the functions take the formatted text), so `GetMaxPrintfLength`
//! (an upper-bound estimator for the C `vsprintf` buffer) is not translated
//! and the buffer growth uses the exact length instead of the estimate; the
//! reading functions return the line instead of filling a caller buffer.
//! `print_sequence_of_nums_compressing_ranges` lives in `ichiprt2.c` and is
//! ported with that file.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust.
//! Licensed under the MIT license (see LICENSES/InChI-MIT.txt).

use std::io::{Read, Write};

use super::util::lrtrim;

pub const INCHI_ADD_STR_LEN: i32 = 32768;
pub const INCHI_STRBUF_INITIAL_SIZE: i32 = 262144;
pub const INCHI_STRBUF_SIZE_INCREMENT: i32 = 262144;
pub const INCHI_STRBUF_SMALLER_INITIAL_SIZE: i32 = 1024;
pub const INCHI_STRBUF_SMALLER_SIZE_INCREMENT: i32 = 4096;

pub const INCHI_IOS_TYPE_NONE: i32 = 0;
pub const INCHI_IOS_TYPE_STRING: i32 = 1;
pub const INCHI_IOS_TYPE_FILE: i32 = 2;

/// `INCHI_IOS_STRING`: `p_str` is the zero-filled allocation (`None` is a
/// NULL pointer), `n_ptr` the read position (or, for a print buffer, the
/// size increment).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InchiIosString {
    pub p_str: Option<Vec<u8>>,
    pub n_allocated_length: i32,
    pub n_used_length: i32,
    pub n_ptr: i32,
}

/// The `FILE *` of a stream
#[derive(Default)]
pub enum IosFile {
    #[default]
    None,
    Stdin,
    Stdout,
    Stderr,
    Reader(Box<dyn Read + Send>),
    Writer(Box<dyn Write + Send>),
}

impl std::fmt::Debug for IosFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            IosFile::None => "None",
            IosFile::Stdin => "Stdin",
            IosFile::Stdout => "Stdout",
            IosFile::Stderr => "Stderr",
            IosFile::Reader(_) => "Reader",
            IosFile::Writer(_) => "Writer",
        })
    }
}

impl IosFile {
    fn is_std(&self) -> bool {
        matches!(self, IosFile::Stdin | IosFile::Stdout | IosFile::Stderr)
    }

    fn getc(&mut self) -> Option<u8> {
        let mut b = [0u8; 1];
        let n = match self {
            IosFile::Stdin => std::io::stdin().read(&mut b).ok()?,
            IosFile::Reader(r) => r.read(&mut b).ok()?,
            _ => return None,
        };
        (n == 1).then_some(b[0])
    }

    fn write(&mut self, data: &[u8]) -> i32 {
        let ok = match self {
            IosFile::Stdout => std::io::stdout().write_all(data).is_ok(),
            IosFile::Stderr => std::io::stderr().write_all(data).is_ok(),
            IosFile::Writer(w) => w.write_all(data).is_ok(),
            _ => return -1,
        };
        if ok { data.len() as i32 } else { -1 }
    }

    fn flush(&mut self) {
        match self {
            IosFile::Stdout => drop(std::io::stdout().flush()),
            IosFile::Stderr => drop(std::io::stderr().flush()),
            IosFile::Writer(w) => drop(w.flush()),
            _ => {}
        }
    }
}

/// `INCHI_IOSTREAM`
#[derive(Debug, Default)]
pub struct InchiIostream {
    pub s: InchiIosString,
    pub f: IosFile,
    pub type_: i32,
}

/// `inchi_ios_init`
pub fn inchi_ios_init(ios: &mut InchiIostream, io_type: i32, f: IosFile) {
    *ios = InchiIostream::default();
    ios.type_ = if io_type == INCHI_IOS_TYPE_FILE {
        INCHI_IOS_TYPE_FILE
    } else {
        INCHI_IOS_TYPE_STRING
    };
    ios.f = f;
}

/// `inchi_ios_create_copy`: copies the type and the buffer size (the
/// content is not copied); a file cannot be shared, so the copy has none.
pub fn inchi_ios_create_copy(ios0: &InchiIostream) -> InchiIostream {
    let mut ios = InchiIostream {
        type_: ios0.type_,
        ..Default::default()
    };
    if ios.type_ == INCHI_IOS_TYPE_STRING {
        ios.s.p_str = Some(vec![0; ios0.s.n_allocated_length.max(0) as usize]);
        ios.s.n_used_length = ios0.s.n_used_length;
        ios.s.n_ptr = ios0.s.n_ptr;
        // the C code leaves nAllocatedLength 0 here
    }
    ios
}

fn free_buffer(s: &mut InchiIosString) {
    s.p_str = None;
    s.n_used_length = 0;
    s.n_allocated_length = 0;
    s.n_ptr = 0;
}

fn used_bytes(s: &InchiIosString) -> &[u8] {
    match &s.p_str {
        Some(p) => &p[..(s.n_used_length.max(0) as usize).min(p.len())],
        None => &[],
    }
}

/// `inchi_ios_flush`: a string buffer is written to its file, then freed
pub fn inchi_ios_flush(ios: &mut InchiIostream) {
    if ios.type_ == INCHI_IOS_TYPE_STRING {
        if ios.s.p_str.is_some() && ios.s.n_used_length > 0 {
            let text = c_prefix(used_bytes(&ios.s)).to_vec();
            if !matches!(ios.f, IosFile::None) {
                ios.f.write(&text);
                ios.f.flush();
            }
            free_buffer(&mut ios.s);
        }
    } else if ios.type_ == INCHI_IOS_TYPE_FILE {
        ios.f.flush();
    }
}

/// The text up to the first NUL (C `%s` of a buffer)
fn c_prefix(b: &[u8]) -> &[u8] {
    &b[..b.iter().position(|&c| c == 0).unwrap_or(b.len())]
}

/// `inchi_ios_flush2`: as [`inchi_ios_flush`], echoing to `f2` as well
pub fn inchi_ios_flush2(ios: &mut InchiIostream, f2: &mut IosFile) {
    if ios.type_ == INCHI_IOS_TYPE_STRING {
        if ios.s.p_str.is_some() && ios.s.n_used_length > 0 {
            let text = c_prefix(used_bytes(&ios.s)).to_vec();
            if !matches!(ios.f, IosFile::None) {
                ios.f.write(&text);
                ios.f.flush();
            }
            f2.write(&text);
            free_buffer(&mut ios.s);
        }
    } else if ios.type_ == INCHI_IOS_TYPE_FILE {
        if !ios.f.is_std() {
            ios.f.flush();
        }
        if !f2.is_std() {
            f2.flush();
        }
    }
}

/// `inchi_ios_close`: frees the buffer and closes the file
pub fn inchi_ios_close(ios: &mut InchiIostream) {
    free_buffer(&mut ios.s);
    if !ios.f.is_std() {
        ios.f = IosFile::None;
    }
}

/// `inchi_ios_reset`: forgets the buffer and closes the file
pub fn inchi_ios_reset(ios: &mut InchiIostream) {
    inchi_ios_close(ios);
}

/// `inchi_ios_free_str`
pub fn inchi_ios_free_str(ios: &mut InchiIostream) {
    free_buffer(&mut ios.s);
}

/// `inchi_ios_str_getc`
pub fn inchi_ios_str_getc(ios: &mut InchiIostream) -> Option<u8> {
    if ios.type_ == INCHI_IOS_TYPE_STRING {
        if ios.s.n_ptr < ios.s.n_used_length {
            let c = ios.s.p_str.as_ref()?.get(ios.s.n_ptr as usize).copied()?;
            ios.s.n_ptr += 1;
            return Some(c);
        }
        return None;
    }
    if ios.type_ == INCHI_IOS_TYPE_FILE {
        return ios.f.getc();
    }
    None
}

/// Common body of `inchi_ios_str_gets` and `inchi_ios_str_getsTab`
fn str_gets_impl(len: i32, ios: &mut InchiIostream, tab_is_lf: bool) -> Option<Vec<u8>> {
    let len = len - 1;
    if len < 0 {
        return None;
    }
    let mut line = Vec::new();
    let mut last = Some(0u8);
    while (line.len() as i32) < len {
        last = inchi_ios_str_getc(ios);
        let Some(mut c) = last else {
            break;
        };
        if tab_is_lf && c == b'\t' {
            c = b'\n';
        }
        line.push(c);
        if c == b'\n' {
            break;
        }
    }
    if line.is_empty() && last.is_none() {
        return None;
    }
    Some(line)
}

/// `inchi_ios_str_gets`: at most `len - 1` bytes up to and including LF
pub fn inchi_ios_str_gets(len: i32, ios: &mut InchiIostream) -> Option<Vec<u8>> {
    str_gets_impl(len, ios, false)
}

/// `inchi_ios_str_getsTab`: as above, a tab ends the line like LF
pub fn inchi_ios_str_gets_tab(len: i32, ios: &mut InchiIostream) -> Option<Vec<u8>> {
    str_gets_impl(len, ios, true)
}

/// The result of the `inchi_ios_gets*` functions: the trimmed length (-1 at
/// the end of the input), the line and whether it was cut at the buffer size.
#[derive(Debug, Clone, PartialEq)]
pub struct IosLine {
    pub length: i32,
    pub line: Vec<u8>,
    pub too_long: bool,
}

fn gets_loop(len: i32, ios: &mut InchiIostream, tab: bool, skip_empty: bool) -> IosLine {
    loop {
        let got = if tab {
            inchi_ios_str_gets_tab(len - 1, ios)
        } else {
            inchi_ios_str_gets(len - 1, ios)
        };
        let Some(mut line) = got else {
            return IosLine {
                length: -1,
                line: Vec::new(),
                too_long: false,
            };
        };
        // szLine[len - 1] = '\0'
        if len >= 1 {
            line.truncate((len - 1) as usize);
        }
        let too_long = !line.contains(&b'\n') && line.len() as i32 == len - 2;
        let length = lrtrim(&mut line);
        if length != 0 || !skip_empty {
            return IosLine {
                length,
                line,
                too_long,
            };
        }
    }
}

/// `inchi_ios_gets`: skips empty lines, trims white space
pub fn inchi_ios_gets(len: i32, ios: &mut InchiIostream) -> IosLine {
    gets_loop(len, ios, false, true)
}

/// `inchi_ios_getsTab`
pub fn inchi_ios_gets_tab(len: i32, ios: &mut InchiIostream) -> IosLine {
    gets_loop(len, ios, true, true)
}

/// `inchi_ios_getsTab1`: as `getsTab` but does not skip empty lines
pub fn inchi_ios_gets_tab1(len: i32, ios: &mut InchiIostream) -> IosLine {
    gets_loop(len, ios, true, false)
}

/// Makes room for `max_len` more bytes in the stream string like the C
/// printing functions; returns false if impossible.
fn ios_ensure(s: &mut InchiIosString, max_len: i32) {
    if s.n_allocated_length - s.n_used_length <= max_len {
        let add = INCHI_ADD_STR_LEN.max(max_len);
        let mut new_str = vec![0u8; (s.n_allocated_length + add) as usize];
        if let Some(old) = &s.p_str {
            let n = s.n_used_length.max(0) as usize;
            new_str[..n].copy_from_slice(&old[..n]);
        }
        s.p_str = Some(new_str);
        s.n_allocated_length += add;
    }
}

fn ios_append(ios: &mut InchiIostream, text: &[u8]) -> i32 {
    let max_len = text.len() as i32;
    ios_ensure(&mut ios.s, max_len);
    let used = ios.s.n_used_length as usize;
    let p = ios.s.p_str.as_mut().expect("buffer was just allocated");
    p[used..used + text.len()].copy_from_slice(text);
    p[used + text.len()] = 0;
    ios.s.n_used_length += max_len;
    max_len
}

/// `inchi_ios_print`: prints formatted text to the stream (-1 on error)
pub fn inchi_ios_print(ios: Option<&mut InchiIostream>, text: &str) -> i32 {
    let Some(ios) = ios else {
        return -1;
    };
    if ios.type_ == INCHI_IOS_TYPE_STRING {
        return ios_append(ios, text.as_bytes());
    }
    if ios.type_ == INCHI_IOS_TYPE_FILE {
        if matches!(ios.f, IosFile::None) {
            return IosFile::Stdout.write(text.as_bytes());
        }
        return ios.f.write(text.as_bytes());
    }
    0
}

/// `inchi_ios_print_nodisplay`: same as print for string streams (the file
/// variant writes to the file, or to stdout if there is none)
pub fn inchi_ios_print_nodisplay(ios: Option<&mut InchiIostream>, text: &str) -> i32 {
    let Some(ios) = ios else {
        return -1;
    };
    if ios.type_ == INCHI_IOS_TYPE_STRING {
        return ios_append(ios, text.as_bytes());
    }
    if ios.type_ == INCHI_IOS_TYPE_FILE {
        if matches!(ios.f, IosFile::None) {
            IosFile::Stdout.write(text.as_bytes());
        } else {
            ios.f.write(text.as_bytes());
        }
    }
    0
}

/// `inchi_ios_flush_not_displayed`: re-prints the buffer content from the
/// start; returns the number of characters, -1 on error
pub fn inchi_ios_flush_not_displayed(ios: Option<&mut InchiIostream>) -> i32 {
    let Some(ios) = ios else {
        return -1;
    };
    let Some(p) = &ios.s.p_str else {
        // C: strcpy from a NULL buffer
        return -1;
    };
    let obuf = c_prefix(&p[..(ios.s.n_used_length.max(0) as usize).min(p.len())]).to_vec();
    ios.s.n_used_length = 0;
    inchi_ios_print(Some(ios), &String::from_utf8_lossy(&obuf))
}

/// `inchi_ios_eprint`: print to the string buffer, or to the file (the
/// stderr echo of the standalone program is not part of the library build)
pub fn inchi_ios_eprint(ios: Option<&mut InchiIostream>, text: &str) -> i32 {
    let Some(ios) = ios else {
        return -1;
    };
    if ios.type_ == INCHI_IOS_TYPE_STRING {
        return ios_append(ios, text.as_bytes());
    }
    if ios.type_ == INCHI_IOS_TYPE_FILE && !matches!(ios.f, IosFile::None) && !text.is_empty() {
        return ios.f.write(text.as_bytes());
    }
    0
}

/// `inchi_sgets`: at most `n - 1` bytes from the string stream (`fgets`
/// emulation, stops at LF and at a zero byte)
pub fn inchi_sgets(n: i32, ios: &mut InchiIostream) -> Option<Vec<u8>> {
    if n <= 0 {
        return None;
    }
    let p_str = ios.s.p_str.as_ref()?;
    let mut inp = ios.s.n_ptr as usize;
    let mut n = n;
    let mut out = Vec::new();
    let mut c = 0u8;
    let mut consumed = 0;
    loop {
        n -= 1;
        if n <= 0 {
            break;
        }
        c = p_str.get(inp).copied().unwrap_or(0);
        inp += 1;
        if c == 0 {
            break;
        }
        consumed += 1;
        out.push(c);
        if c == b'\n' {
            break;
        }
    }
    ios.s.n_ptr += consumed;
    if c == 0 && out.is_empty() {
        None
    } else {
        Some(out)
    }
}

/// `inchi_fgetsLf`: reads up to LF but at most `line_len - 1` bytes; the rest
/// of a longer line is skipped. A CR is turned into LF and ends the line.
pub fn inchi_fgets_lf(line_len: i32, inp_stream: &mut InchiIostream) -> Option<Vec<u8>> {
    let mut p: Option<Vec<u8>> = None;
    if inp_stream.type_ == INCHI_IOS_TYPE_FILE {
        // fgets
        if line_len > 1 {
            let mut buf = Vec::new();
            while (buf.len() as i32) < line_len - 1 {
                match inp_stream.f.getc() {
                    Some(c) => {
                        buf.push(c);
                        if c == b'\n' {
                            break;
                        }
                    }
                    None => break,
                }
            }
            if !buf.is_empty() {
                p = Some(buf);
            }
        } else if line_len == 1 {
            p = Some(Vec::new());
        }
        if let Some(l) = &p
            && !l.contains(&b'\n')
        {
            // bypass up to LF or the end of the file
            while let Some(c) = inp_stream.f.getc() {
                if c == b'\n' {
                    break;
                }
            }
        }
    } else if inp_stream.type_ == INCHI_IOS_TYPE_STRING {
        p = inchi_sgets(line_len, inp_stream);
        if let Some(l) = &p
            && !l.contains(&b'\n')
        {
            while let Some(t) = inchi_sgets(64, inp_stream) {
                if t.contains(&b'\n') {
                    break;
                }
            }
        }
    }
    if let Some(l) = &mut p
        && let Some(q) = l.iter().position(|&c| c == b'\r')
    {
        // fix CR CR LF line terminator
        l.truncate(q);
        l.push(b'\n');
    }
    p
}

// ---- the growing print buffer ----

/// `inchi_strbuf_init`: returns the start size
pub fn inchi_strbuf_init(buf: &mut InchiIosString, start_size: i32, incr_size: i32) -> i32 {
    *buf = InchiIosString::default();
    let start_size = if start_size <= 0 {
        INCHI_STRBUF_INITIAL_SIZE
    } else {
        start_size
    };
    let incr_size = if incr_size <= 0 {
        INCHI_STRBUF_SIZE_INCREMENT
    } else {
        incr_size
    };
    buf.p_str = Some(vec![0; start_size as usize]);
    buf.n_allocated_length = start_size;
    buf.n_ptr = incr_size;
    start_size
}

/// `inchi_strbuf_reset`: empties the string; note that the increment (kept in
/// `n_ptr`) is reset too, as in the C code
pub fn inchi_strbuf_reset(buf: &mut InchiIosString) {
    if let Some(p) = &mut buf.p_str
        && !p.is_empty()
    {
        p[0] = 0;
    }
    buf.n_used_length = 0;
    buf.n_ptr = 0;
}

/// `inchi_strbuf_close`
pub fn inchi_strbuf_close(buf: &mut InchiIosString) {
    *buf = InchiIosString::default();
}

/// `inchi_strbuf_create_copy`: same sizes, zeroed content
pub fn inchi_strbuf_create_copy(buf: &InchiIosString) -> InchiIosString {
    InchiIosString {
        p_str: Some(vec![0; buf.n_allocated_length.max(0) as usize]),
        n_allocated_length: buf.n_allocated_length,
        n_used_length: buf.n_used_length,
        n_ptr: buf.n_ptr,
    }
}

/// `inchi_strbuf_update`: makes room for `new_addition_size` more bytes;
/// returns the allocated length
pub fn inchi_strbuf_update(buf: &mut InchiIosString, new_addition_size: i32) -> i32 {
    if new_addition_size <= 0 {
        return buf.n_allocated_length;
    }
    let requested = buf.n_used_length + new_addition_size;
    if requested >= buf.n_allocated_length {
        let add = buf.n_ptr.max(new_addition_size);
        let mut new_str = vec![0u8; (buf.n_allocated_length + add) as usize];
        if let Some(old) = &buf.p_str {
            let n = buf.n_used_length.max(0) as usize;
            new_str[..n].copy_from_slice(&old[..n]);
        }
        buf.p_str = Some(new_str);
        buf.n_allocated_length += add;
    }
    buf.n_allocated_length
}

/// `inchi_strbuf_printf`: appends the text; returns its length
pub fn inchi_strbuf_printf(buf: &mut InchiIosString, text: &str) -> i32 {
    let b = text.as_bytes();
    inchi_strbuf_update(buf, b.len() as i32);
    let used = buf.n_used_length as usize;
    let p = buf.p_str.as_mut().expect("buffer was allocated by update");
    p[used..used + b.len()].copy_from_slice(b);
    p[used + b.len()] = 0;
    buf.n_used_length += b.len() as i32;
    b.len() as i32
}

/// `inchi_strbuf_printf_from`: writes the text at position `npos` (no
/// checks, as in C) and sets the used length to `npos + len`
pub fn inchi_strbuf_printf_from(buf: &mut InchiIosString, npos: i32, text: &str) -> i32 {
    let b = text.as_bytes();
    inchi_strbuf_update(buf, b.len() as i32 + npos);
    let pos = npos as usize;
    let p = buf.p_str.as_mut().expect("buffer was allocated by update");
    p[pos..pos + b.len()].copy_from_slice(b);
    p[pos + b.len()] = 0;
    buf.n_used_length = npos + b.len() as i32;
    b.len() as i32
}

/// `inchi_strbuf_printf(buf, "%c", c)`: `GetMaxPrintfLength` estimates two
/// bytes for `%c`, which decides when the buffer grows.
fn strbuf_put_char(buf: &mut InchiIosString, c: u8) {
    inchi_strbuf_update(buf, 2);
    let used = buf.n_used_length as usize;
    let p = buf.p_str.as_mut().expect("buffer was allocated by update");
    p[used] = c;
    p[used + 1] = 0;
    buf.n_used_length += 1;
}

fn finish_line(buf: &mut InchiIosString, crlf2lf: bool, preserve_lf: bool) -> i32 {
    let p = buf.p_str.as_mut().expect("line was added");
    if crlf2lf && buf.n_used_length > 2 && p[(buf.n_used_length - 2) as usize] == b'\r' {
        p[(buf.n_used_length - 2) as usize] = b'\n';
        buf.n_used_length -= 1;
        p[buf.n_used_length as usize] = 0;
    }
    if !preserve_lf {
        buf.n_used_length -= 1;
        p[buf.n_used_length as usize] = 0;
    }
    buf.n_used_length
}

/// `inchi_strbuf_getline`: resets the buffer and reads the next line of the
/// file into it; the number of characters, -1 at the end or on error
pub fn inchi_strbuf_getline(
    buf: &mut InchiIosString,
    f: &mut IosFile,
    crlf2lf: bool,
    preserve_lf: bool,
) -> i32 {
    inchi_strbuf_reset(buf);
    loop {
        let Some(c) = f.getc() else {
            return -1;
        };
        strbuf_put_char(buf, c);
        if c == b'\n' {
            break;
        }
    }
    finish_line(buf, crlf2lf, preserve_lf)
}

/// `inchi_strbuf_addline`: adds the next line of the stream to the buffer
pub fn inchi_strbuf_addline(
    buf: &mut InchiIosString,
    inp_stream: &mut InchiIostream,
    crlf2lf: bool,
    preserve_lf: bool,
) -> i32 {
    loop {
        let Some(c) = inchi_ios_str_getc(inp_stream) else {
            return -1;
        };
        strbuf_put_char(buf, c);
        if c == b'\n' {
            break;
        }
    }
    finish_line(buf, crlf2lf, preserve_lf)
}

/// `Output_RecordInfo`: the header of the current record
pub fn output_record_info(
    out_file: &mut InchiIostream,
    num_input_struct: i32,
    b_no_struct_labels: bool,
    sz_sdf_label: Option<&str>,
    sz_sdf_value: Option<&str>,
    l_sdf_id: u64,
    p_lf: &str,
    p_tab: &str,
) -> i32 {
    if b_no_struct_labels {
        return 0;
    }
    let label = sz_sdf_label.filter(|s| !s.is_empty());
    let value = sz_sdf_value.filter(|s| !s.is_empty());
    if label.is_none() && value.is_none() {
        inchi_ios_print_nodisplay(
            Some(out_file),
            &format!("{p_lf}Structure: {num_input_struct}"),
        );
        inchi_ios_print_nodisplay(Some(out_file), p_tab);
    } else {
        // SDF_LBL_VAL
        let (a, b, c, d) = match (label, value) {
            (Some(l), Some(v)) => (" ", l, "=", v),
            (Some(l), None) => (" ", l, " ", "is missing"),
            (None, v) => ("", "", "", v.unwrap_or("")),
        };
        inchi_ios_print_nodisplay(
            Some(out_file),
            &format!("{p_lf}Structure: {num_input_struct}.{a}{b}{c}{d}"),
        );
        if l_sdf_id != 0 {
            out_file.s.n_used_length -= 1;
            inchi_ios_print_nodisplay(Some(out_file), &format!(":{l_sdf_id}"));
        }
        inchi_ios_print_nodisplay(Some(out_file), p_tab);
    }
    0
}
