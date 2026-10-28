# Plan: PNG drawing (`MolDraw2DCairo`) without Cairo

## Where we are

The whole `MolDraw2D` stack is ported (stages 53–71): `DrawMol`, shapes, text (FreeType glyph outlines read by our own TrueType
reader), highlights, reactions, contours, and the SVG drawer (`MolDraw2DSVG`). Only the raster back end is missing:
`MolDraw2DCairo`, `Draw.MolToImage`/`MolsToGridImage` as PNG, `MolDraw2D.GetDrawingText()` bytes, and the upstream tests
that need PNG (`UnitTestDraw.py` 11/16, parts of `testMolDraw2D.py`). PNG *metadata* read/write already exists
(`file_parsers::png`, stage 120).

## Can Cairo itself be used?

`cairo-rs` is only FFI to the C library `libcairo` (plus pixman, freetype, fontconfig). It needs a system library or a bundled C
build, breaks the pure-Rust/`unsafe_code = forbid`/abi3-wheel approach of this repository and makes wheels platform-fragile.
Not recommended. Pixel-identical output with Cairo would need porting cairo's stroker and pixman's rasteriser (tens of
thousands of lines) and is not worth it.

## Ready-made pure Rust options

| crate | what | verdict |
|---|---|---|
| `tiny-skia` (BSD-3, pure Rust, no unsafe needed on our side) | anti-aliased fill/stroke of paths, caps, joins, dashes, clip, alpha; PNG encode via `png` | **recommended**: same drawing model as Cairo (path, fill, stroke) |
| `png` + `miniz_oxide` | encoder | pulled in by tiny-skia; our own metadata writer stays for tEXt/iTXt chunks |
| `resvg` (renders our SVG to PNG) | cheapest, but heavy dependency tree (usvg, fontdb, ...) and no control of the drawing calls | alternative if a quick path is wanted |
| `raqote`, `vello_cpu`, `plotters` | similar, less mature or bigger | not chosen |

Text needs no font library: the drawer already turns text into glyph outlines (`DrawTextFT` paths) from the built-in Telex/
Roboto fonts; the raster back end just fills those paths. `noFreetype=True` (Cairo toy fonts) falls back to the same outlines.

## Stages (one branch each, tests, PR, merge when CI is green)

1. **P1 raster canvas** (`mol_draw2d::raster_canvas`): a `ShapeCanvas`/drawer back end over `tiny-skia::Pixmap`: background,
   lines (butt caps, scaled width, dashes), polygons (fill/outline), ellipses, arcs, wavy lines, glyph-path text, alpha.
   Behind a cargo feature `png` (on in the Python wheels) so the core crate keeps zero drawing dependencies.
2. **P2 `MolDraw2DCairo` API**: constructor(w, h, panelW, panelH, noFreetype), `drawMolecule(s)`, `finishDrawing`,
   `getDrawingText()` (PNG bytes), `writeDrawingText(path)`, `addMoleculeMetadata` (reuses `file_parsers::png`),
   `MolToPNG`-style helpers, reaction drawing.
3. **P3 Python**: `rdMolDraw2D.MolDraw2DCairo`, `Draw.MolToImage`, `MolToFile`, `MolsToGridImage(useSVG=False)`,
   `DrawMorganBit`/`DrawRDKitBit` PNG, IPython PNG display. PIL stays an optional dependency (pure Python upstream code).
4. **P4 verification and ratchet**: upstream `UnitTestDraw.py` PNG cases and the Cairo cases of `testMolDraw2D.py`; direct
   comparison with the official wheel's PNG using image metrics.
5. **P5 (optional)**: heat maps/`SimilarityMaps` PNG once matplotlib is not required; JPEG is out of scope.

## What "identical" means here

Geometry, colours, sizes, text layout, API and PNG metadata can be identical to RDKit. Anti-aliased edge pixels cannot be
bit-identical to Cairo/pixman. Tests therefore check: image size, ink bounding box, exact colours at known interior pixels,
metadata round trip, and a tolerance on mean absolute difference (and SSIM) against the official wheel on a corpus; the
tolerance is recorded in the test, not hidden.

## Size and risks

About 2–3k lines of Rust for P1–P2, small for P3. Risks: tiny-skia stroke joins/caps differ slightly from Cairo's
(mitigated by drawing bonds as filled polygons where RDKit already does); dependency growth (feature-gated, `cargo deny`
license check: tiny-skia/png/miniz_oxide are BSD/MIT/Apache); wheel size (+~1 MB).
