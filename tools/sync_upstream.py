#!/usr/bin/env python3
"""Keep rdoxide in sync with an upstream RDKit release.

Usage (Python >= 3.11, stdlib only):

    python3 tools/sync_upstream.py fetch    # clone/checkout UPSTREAM.toml tag into upstream/rdkit
    python3 tools/sync_upstream.py regen    # re-extract data tables and copy pure-python files/test data
    python3 tools/sync_upstream.py diff     # upstream files changed since they were last ported
    python3 tools/sync_upstream.py stamp [MODULE...]  # record current upstream blob hashes as ported

`porting/manifest.toml` maps every Rust module to the upstream files it was
ported from, together with the git blob hash of each file at port time. When a
new RDKit release is checked out, `diff` lists exactly which ported upstream
files changed so the corresponding Rust code can be updated.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
UPSTREAM_DIR = ROOT / "upstream" / "rdkit"
MANIFEST = ROOT / "porting" / "manifest.toml"
DATA = ROOT / "crates" / "rdoxide" / "data"
COPIES = ROOT / "porting" / "copies.toml"


def load_upstream_cfg() -> dict:
    with open(ROOT / "UPSTREAM.toml", "rb") as fh:
        return tomllib.load(fh)["rdkit"]


def git(*args: str, cwd: Path = UPSTREAM_DIR) -> str:
    return subprocess.run(
        ["git", *args], cwd=cwd, check=True, capture_output=True, text=True
    ).stdout.strip()


# --------------------------------------------------------------------------- fetch
def cmd_fetch(_args) -> None:
    cfg = load_upstream_cfg()
    if not UPSTREAM_DIR.exists():
        UPSTREAM_DIR.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["git", "clone", "--depth", "1", "--branch", cfg["tag"], cfg["repository"],
             str(UPSTREAM_DIR)],
            check=True,
        )
    else:
        git("fetch", "--depth", "1", "origin", "tag", cfg["tag"])
        git("checkout", "--quiet", cfg["tag"])
    head = git("rev-parse", "HEAD")
    if head != cfg["commit"]:
        sys.exit(f"upstream HEAD {head} does not match UPSTREAM.toml commit {cfg['commit']}")
    print(f"upstream/rdkit at {cfg['tag']} ({head})")


# --------------------------------------------------------------------------- regen
_RAW_RE = re.compile(r'R"DAT\((.*?)\)DAT"', re.S)


def _cxx_raw_concat(src: str, start_marker: str) -> str:
    """Concatenate the R"DAT(...)DAT" literals of the initializer after start_marker."""
    # drop whole-line // comments (they may contain ';')
    src = re.sub(r"(?m)^\s*//.*$", "", src)
    start = pos = src.index(start_marker)
    while True:
        # the initializer ends at the first ';' that is not inside a raw literal
        lit = src.find('R"DAT(', pos)
        semi = src.find(";", pos)
        if lit == -1 or semi < lit:
            return "".join(_RAW_RE.findall(src[start:semi]))
        pos = src.index(')DAT"', lit) + 5


def regen_atomic_data() -> None:
    src = (UPSTREAM_DIR / "Code/GraphMol/atomic_data.cpp").read_text()
    header = (
        "# Extracted from RDKit Code/GraphMol/atomic_data.cpp by tools/sync_upstream.py.\n"
        "# Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors (BSD-3-Clause).\n"
        "# DO NOT EDIT - regenerate with `python3 tools/sync_upstream.py regen`.\n"
    )
    # periodic table: one element per line
    table = _cxx_raw_concat(src, "periodicTableAtomData =")
    lines = [ln.strip() for ln in table.splitlines() if ln.strip()]
    (DATA / "atomic_data.tsv").write_text(
        header
        + "# anum symbol row rCov rB0 rVdw mass nOuterElecs commonIsotope commonIsotopeMass valences...\n"
        + "\n".join("\t".join(ln.split()) for ln in lines)
        + "\n"
    )
    # element names
    m = re.search(r"elementNames\s*=\s*\{(.*?)\};", src, re.S)
    names = re.findall(r'"([^"]*)"', m.group(1))
    (DATA / "element_names.txt").write_text(header + "\n".join(names) + "\n")
    # isotopes (lines in upstream end with a literal "\n \" sequence)
    iso = _cxx_raw_concat(src, "isotopesAtomData[] =")
    out = []
    for ln in iso.splitlines():
        ln = ln.replace("\\n", " ").replace("\\", " ").strip()
        if ln:
            out.append("\t".join(ln.split()[:5]))
    (DATA / "isotopes.tsv").write_text(
        header + "# anum symbol isotope mass abundance\n" + "\n".join(out) + "\n"
    )
    print(f"atomic_data.tsv: {len(lines)} rows; data/isotopes.tsv: {len(out)} rows; "
          f"data/element_names.txt: {len(names)} names")


def regen_common_properties() -> None:
    """Generate rdgeneral/common_properties.rs from Code/RDGeneral/types.h."""
    src = (UPSTREAM_DIR / "Code/RDGeneral/types.h").read_text()
    block = src[src.index("namespace common_properties {"):]
    block = block[: block.index("}  // namespace common_properties")]
    out = [
        "//! Property names shared across RDKit (`RDKit::common_properties`).",
        "//!",
        "//! GENERATED from RDKit `Code/RDGeneral/types.h` by tools/sync_upstream.py regen.",
        "//! DO NOT EDIT. Copyright (C) RDKit contributors, BSD-3-Clause.",
        "#![allow(non_upper_case_globals)]",
        "",
    ]
    # declarations may wrap onto a second line
    block = re.sub(r"//[^\n]*", "", block)
    tokens = re.findall(
        r'inline constexpr std::string_view (\w+)\s*=\s*("(?:[^"\\]|\\.)*");'
        r"|namespace (\w+) \{|(\})",
        block,
    )
    n, depth = 0, 0
    for name, value, ns, close in tokens:
        if name:
            out.append(f"{'    ' * depth}pub const {name}: &str = {value};")
            n += 1
        elif ns == "common_properties":
            continue
        elif ns:
            out.append(f"{'    ' * depth}pub mod {ns} {{")
            depth += 1
        elif close and depth:
            depth -= 1
            out.append(f"{'    ' * depth}}}")
    path = ROOT / "crates/rdoxide/src/rdgeneral/common_properties.rs"
    path.write_text("\n".join(out) + "\n")
    print(f"common_properties.rs: {n} names")


def regen_copies() -> None:
    """Copy verbatim files (pure-python modules, test data) listed in porting/copies.toml."""
    if not COPIES.exists():
        return
    with open(COPIES, "rb") as fh:
        cfg = tomllib.load(fh)
    n = 0
    for entry in cfg.get("copy", []):
        src = UPSTREAM_DIR / entry["from"]
        dst = ROOT / entry["to"]
        if src.is_dir():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst, ignore=shutil.ignore_patterns(*entry.get("exclude", [])))
        else:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
        n += 1
    print(f"copied {n} upstream paths")


def cmd_regen(_args) -> None:
    if not UPSTREAM_DIR.exists():
        sys.exit("upstream/rdkit missing - run `sync_upstream.py fetch` first")
    regen_atomic_data()
    regen_common_properties()
    regen_copies()


# --------------------------------------------------------------------------- diff / stamp
def _load_manifest() -> dict:
    with open(MANIFEST, "rb") as fh:
        return tomllib.load(fh)


def _blob(path: str) -> str:
    return git("rev-parse", f"HEAD:{path}")


def cmd_diff(_args) -> None:
    manifest = _load_manifest()
    stale = 0
    for module, spec in manifest.get("module", {}).items():
        for path, blob in spec.get("sources", {}).items():
            try:
                cur = _blob(path)
            except subprocess.CalledProcessError:
                cur = "<deleted>"
            if cur != blob:
                stale += 1
                print(f"{module}: {path} changed ({blob[:10]} -> {cur[:10]})")
    print(f"{stale} upstream file(s) changed since last port" if stale else "all ported files up to date")
    sys.exit(1 if stale else 0)


def cmd_stamp(args) -> None:
    """Rewrite blob hashes in manifest.toml (simple line-based rewrite keeps comments)."""
    text = MANIFEST.read_text()
    out, module = [], None
    for line in text.splitlines():
        mh = re.match(r'\[module\."([^"]+)"', line) or re.match(r"\[module\.([\w.]+)", line)
        if mh:
            module = mh.group(1)
        ms = re.match(r'^(\s*)"([^"]+)"\s*=\s*"[^"]*"(.*)$', line)
        if ms and module and (not args.modules or module in args.modules):
            indent, path, rest = ms.groups()
            line = f'{indent}"{path}" = "{_blob(path)}"{rest}'
        out.append(line)
    MANIFEST.write_text("\n".join(out) + "\n")
    print("manifest stamped")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("fetch").set_defaults(fn=cmd_fetch)
    sub.add_parser("regen").set_defaults(fn=cmd_regen)
    sub.add_parser("diff").set_defaults(fn=cmd_diff)
    s = sub.add_parser("stamp")
    s.add_argument("modules", nargs="*")
    s.set_defaults(fn=cmd_stamp)
    args = p.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
