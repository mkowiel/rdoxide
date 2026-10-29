"""Drawing cases shared by tools/make_png_reference.py (official wheel) and test_png_reference.py (rdoxide)."""

CASES = {
    "aspirin": ("CC(=O)Oc1ccccc1C(=O)O", (300, 250), {}),
    "caffeine": ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", (300, 300), {}),
    "stereo": ("C[C@H](N)C(=O)O", (250, 250), {}),
    "ez": ("C/C=C/C(=O)O", (250, 200), {}),
    "charged": ("C[N+](C)(C)CC(=O)[O-]", (250, 250), {}),
    "steroid": ("CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O", (400, 300), {}),
    "legend": ("c1ccncc1", (250, 250), {"legend": "pyridine"}),
    "highlight": ("c1ccccc1CCO", (300, 250), {"highlight": [0, 1, 2]}),
    "macro": ("C1COCCOCCOCCOCCOCCO1", (300, 300), {}),
    "bigger": ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", (500, 400), {}),
}


def draw(Chem, rdMolDraw2D, name):
    smi, (w, h), opt = CASES[name]
    m = Chem.MolFromSmiles(smi)
    d = rdMolDraw2D.MolDraw2DCairo(w, h)
    if "highlight" in opt:
        d.DrawMolecule(m, highlightAtoms=opt["highlight"])
    else:
        d.DrawMolecule(m, legend=opt.get("legend", ""))
    d.FinishDrawing()
    return d.GetDrawingText()
