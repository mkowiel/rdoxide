//! Ported from RDKit `Code/RDGeneral/RDLog.{h,cpp}`.
//!
//! Copyright (C) 2005-2022 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! Four global logs (`rdApp.debug`, `rdApp.info`, `rdApp.warning`,
//! `rdApp.error`). Every message is prefixed with a `[HH:MM:SS] ` local-time
//! stamp and terminated by a newline, exactly like upstream's `BOOST_LOG`.
//! Destinations are pluggable so the Python layer can route logs to Python's
//! `logging` module or `sys.stderr`.

use std::io::Write;
use std::sync::{Arc, Mutex, MutexGuard};

/// Log channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Level {
    Debug = 0,
    Info = 1,
    Warning = 2,
    Error = 3,
}

impl Level {
    pub const ALL: [Level; 4] = [Level::Debug, Level::Info, Level::Warning, Level::Error];

    pub fn spec(self) -> &'static str {
        match self {
            Level::Debug => "rdApp.debug",
            Level::Info => "rdApp.info",
            Level::Warning => "rdApp.warning",
            Level::Error => "rdApp.error",
        }
    }

    pub fn from_spec(spec: &str) -> Option<Level> {
        Level::ALL.into_iter().find(|l| l.spec() == spec)
    }
}

/// A log destination. `text` is the complete formatted message including the
/// timestamp and trailing newline.
pub trait LogSink: Send + Sync {
    fn write(&self, level: Level, text: &str);
}

/// Writes to the process' stderr (`std::cerr`).
#[derive(Debug)]
pub struct StderrSink;
impl LogSink for StderrSink {
    fn write(&self, _level: Level, text: &str) {
        let _ = std::io::stderr().lock().write_all(text.as_bytes());
    }
}

/// Writes to the process' stdout (`std::cout`).
#[derive(Debug)]
pub struct StdoutSink;
impl LogSink for StdoutSink {
    fn write(&self, _level: Level, text: &str) {
        let mut out = std::io::stdout().lock();
        let _ = out.write_all(text.as_bytes());
        let _ = out.flush();
    }
}

struct Logger {
    enabled: bool,
    dest: Option<Arc<dyn LogSink>>,
    tee: Option<Arc<dyn LogSink>>,
    /// Stack of capture buffers (`CaptureLog`); innermost capture wins.
    captures: Vec<(u64, String)>,
}

impl Logger {
    const fn uninit() -> Self {
        Logger {
            enabled: false,
            dest: None,
            tee: None,
            captures: Vec::new(),
        }
    }
}

struct LogState {
    initialized: bool,
    logs: [Logger; 4],
    next_capture_id: u64,
}

static STATE: Mutex<LogState> = Mutex::new(LogState {
    initialized: false,
    logs: [
        Logger::uninit(),
        Logger::uninit(),
        Logger::uninit(),
        Logger::uninit(),
    ],
    next_capture_id: 1,
});

fn state() -> MutexGuard<'static, LogState> {
    let mut st = STATE.lock().unwrap_or_else(|e| e.into_inner());
    if !st.initialized {
        init_logs_locked(&mut st);
    }
    st
}

fn init_logs_locked(st: &mut LogState) {
    st.initialized = true;
    let stderr: Arc<dyn LogSink> = Arc::new(StderrSink);
    let stdout: Arc<dyn LogSink> = Arc::new(StdoutSink);
    for (lvl, (dest, enabled)) in [
        (Level::Debug, (stderr.clone(), false)),
        (Level::Info, (stdout, false)),
        (Level::Warning, (stderr.clone(), true)),
        (Level::Error, (stderr, true)),
    ] {
        let log = &mut st.logs[lvl as usize];
        log.dest = Some(dest);
        log.enabled = enabled;
        log.tee = None;
    }
}

/// `RDLog::InitLogs()` (Python: `rdBase.LogToCppStreams()`).
pub fn init_logs() {
    let mut st = STATE.lock().unwrap_or_else(|e| e.into_inner());
    init_logs_locked(&mut st);
}

/// Replace the destination of every log, keeping enabled flags.
pub fn set_destinations(make: impl Fn(Level) -> Arc<dyn LogSink>) {
    let mut st = state();
    for lvl in Level::ALL {
        let log = &mut st.logs[lvl as usize];
        log.dest = Some(make(lvl));
        log.tee = None;
    }
}

/// Tee every log to an additional sink (`rdLogger::SetTee`).
pub fn set_tees(make: impl Fn(Level) -> Arc<dyn LogSink>) {
    let mut st = state();
    for lvl in Level::ALL {
        st.logs[lvl as usize].tee = Some(make(lvl));
    }
}

fn levels_for_spec(spec: &str) -> Vec<Level> {
    if spec == "rdApp.*" {
        Level::ALL.to_vec()
    } else {
        Level::from_spec(spec).into_iter().collect()
    }
}

/// `boost::logging::enable_logs`
pub fn enable_logs(spec: &str) {
    let mut st = state();
    for l in levels_for_spec(spec) {
        st.logs[l as usize].enabled = true;
    }
}

/// `boost::logging::disable_logs`
pub fn disable_logs(spec: &str) {
    let mut st = state();
    for l in levels_for_spec(spec) {
        st.logs[l as usize].enabled = false;
    }
}

pub fn is_enabled(level: Level) -> bool {
    state().logs[level as usize].enabled
}

/// `boost::logging::log_status`
pub fn log_status() -> String {
    let st = state();
    Level::ALL
        .iter()
        .map(|&l| {
            let log = &st.logs[l as usize];
            let status = if log.dest.is_none() {
                "unitialized"
            } else if log.enabled {
                "enabled"
            } else {
                "disabled"
            };
            format!("{}:{status}", l.spec())
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn timestamp() -> String {
    chrono::Local::now().format("[%H:%M:%S] ").to_string()
}

/// Emit `msg` (a newline is appended) on `level` if that log is enabled.
pub fn log(level: Level, msg: &str) {
    let (sinks, captured) = {
        let mut st = state();
        let lg = &mut st.logs[level as usize];
        if !lg.enabled || lg.dest.is_none() {
            return;
        }
        let text = format!("{}{msg}\n", timestamp());
        if let Some((_, buf)) = lg.captures.last_mut() {
            buf.push_str(&text);
            (Vec::new(), true)
        } else {
            let mut sinks: Vec<Arc<dyn LogSink>> = Vec::with_capacity(2);
            sinks.extend(lg.dest.clone());
            sinks.extend(lg.tee.clone());
            (
                sinks
                    .into_iter()
                    .map(|s| (s, text.clone()))
                    .collect::<Vec<_>>(),
                false,
            )
        }
    };
    if captured {
        return;
    }
    // call sinks without holding the lock: Python sinks may re-enter
    for (sink, text) in sinks {
        sink.write(level, &text);
    }
}

#[macro_export]
macro_rules! rd_log_error { ($($t:tt)*) => { $crate::rdgeneral::log::log($crate::rdgeneral::log::Level::Error, &format!($($t)*)) } }
#[macro_export]
macro_rules! rd_log_warning { ($($t:tt)*) => { $crate::rdgeneral::log::log($crate::rdgeneral::log::Level::Warning, &format!($($t)*)) } }
#[macro_export]
macro_rules! rd_log_info { ($($t:tt)*) => { $crate::rdgeneral::log::log($crate::rdgeneral::log::Level::Info, &format!($($t)*)) } }
#[macro_export]
macro_rules! rd_log_debug { ($($t:tt)*) => { $crate::rdgeneral::log::log($crate::rdgeneral::log::Level::Debug, &format!($($t)*)) } }

/// Disables all logs while alive and restores the previous state on drop
/// (`RDLog::LogStateSetter`, Python `rdBase.BlockLogs`).
#[derive(Debug)]
pub struct LogStateSetter {
    orig_state: u8,
}

impl LogStateSetter {
    pub fn new() -> Self {
        Self::with_enabled(&[])
    }

    /// Disable all logs except `to_enable`, which are switched on.
    pub fn with_enabled(to_enable: &[Level]) -> Self {
        let mut st = state();
        let mut orig = 0u8;
        for l in Level::ALL {
            let lg = &mut st.logs[l as usize];
            if lg.enabled {
                orig |= 1 << l as u8;
                lg.enabled = false;
            }
        }
        for &l in to_enable {
            orig ^= 1 << l as u8;
            st.logs[l as usize].enabled = true;
        }
        LogStateSetter { orig_state: orig }
    }
}

impl Default for LogStateSetter {
    fn default() -> Self {
        Self::new()
    }
}

impl Drop for LogStateSetter {
    fn drop(&mut self) {
        let mut st = state();
        for l in Level::ALL {
            st.logs[l as usize].enabled ^= (self.orig_state >> l as u8) & 1 == 1;
        }
    }
}

/// Captures messages of one log while alive (`RDLog::CaptureLog`).
#[derive(Debug)]
pub struct CaptureLog {
    level: Level,
    id: u64,
    was_enabled: bool,
    finished: Option<String>,
}

impl CaptureLog {
    pub fn new(level: Level) -> Self {
        let mut st = state();
        let id = st.next_capture_id;
        st.next_capture_id += 1;
        let lg = &mut st.logs[level as usize];
        let was_enabled = lg.enabled;
        lg.enabled = true;
        lg.captures.push((id, String::new()));
        CaptureLog {
            level,
            id,
            was_enabled,
            finished: None,
        }
    }

    /// `RDLog::CaptureErrorLog`
    pub fn error() -> Self {
        Self::new(Level::Error)
    }

    pub fn messages(&self) -> String {
        if let Some(m) = &self.finished {
            return m.clone();
        }
        let st = state();
        st.logs[self.level as usize]
            .captures
            .iter()
            .find(|(id, _)| *id == self.id)
            .map(|(_, m)| m.clone())
            .unwrap_or_default()
    }

    /// Stop capturing; messages remain available.
    pub fn finish(&mut self) {
        if self.finished.is_some() {
            return;
        }
        let mut st = state();
        let lg = &mut st.logs[self.level as usize];
        let msgs = match lg.captures.iter().position(|(id, _)| *id == self.id) {
            Some(pos) => lg.captures.remove(pos).1,
            None => String::new(),
        };
        lg.enabled = self.was_enabled;
        self.finished = Some(msgs);
    }
}

impl Drop for CaptureLog {
    fn drop(&mut self) {
        self.finish();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // the log state is global, so everything is exercised in one test
    #[test]
    fn enable_block_capture() {
        disable_logs("rdApp.*");
        assert_eq!(
            log_status(),
            "rdApp.debug:disabled\nrdApp.info:disabled\nrdApp.warning:disabled\nrdApp.error:disabled"
        );
        enable_logs("rdApp.error");
        assert!(is_enabled(Level::Error));
        {
            let _b = LogStateSetter::new();
            assert!(!is_enabled(Level::Error));
        }
        assert!(is_enabled(Level::Error));
        assert!(!is_enabled(Level::Warning));

        let mut cap = CaptureLog::error();
        log(Level::Error, "boom");
        {
            let inner = CaptureLog::error();
            log(Level::Error, "inner");
            assert!(inner.messages().ends_with("] inner\n"));
        }
        cap.finish();
        let m = cap.messages();
        assert!(m.starts_with('[') && m.ends_with("] boom\n"), "{m}");
        enable_logs("rdApp.*");
    }
}
