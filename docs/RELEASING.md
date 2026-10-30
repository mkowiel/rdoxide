# Releasing

Versions are `<RDKit version>[pre-release]`: the first alpha of the port of
RDKit 2026.09.1 is `2026.9.1a0` (PEP 440, as pip sees it) = `2026.9.1-alpha.0`
(semver, in `Cargo.toml`; maturin converts it).

1. Make sure `main` is green, then run `tools/regen_speed_report.sh` if the
   benchmarks should be refreshed and `python tools/run_upstream_tests.py --update`.
2. Set `version` and the `rdoxide` path dependency version in the workspace
   `Cargo.toml` (and `cargo update -w` for `Cargo.lock`), add a section to
   `CHANGELOG.md` headed `## <pip version>`, and update the numbers in
   `README.md`.
3. `python tools/check_licenses.py` must pass (`--write` after dependency
   changes); `NOTICE`, `LICENSES/` and the SPDX expression in `Cargo.toml` and
   `pyproject.toml` must list every ported third-party component.
4. Merge, then tag `v<pip version>` on `main` (`git tag v2026.9.1a0 && git push
   origin v2026.9.1a0`). The `Release` workflow builds wheels for Linux
   x86-64/aarch64 (manylinux 2.28), macOS arm64/x86-64 and Windows x86-64, the
   sdist, and creates the GitHub release (marked pre-release for a/b/rc tags)
   with the changelog section as notes.
5. PyPI: download the release files and upload them with
   `twine upload dist/*` (or configure trusted publishing for the repository and
   add a publish job). Pre-releases are only installed with `pip install --pre`.
