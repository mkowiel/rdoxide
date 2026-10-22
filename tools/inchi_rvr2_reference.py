"""Generates tests/data/inchi_rvr2.json from the InChI C code: the structure
restoration of ichirvr2.c from the components of the InChI strings of the
official RDKit corpus (`RestoreAtomConnectionsSetStereo`,
`SetStereoBondTypesFrom0DStereo` and their helpers).

    python tools/inchi_rvr2_reference.py
"""
import importlib.util
import json
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


r1 = load("rvr1", "tools/inchi_rvr1_reference.py")
r2, ma, SRC, hx = r1.r2, r1.ma, r1.SRC, r1.hx

MAIN = r1.PRE + r"""
static void dump_atoms(StrFromINChI *s) {
    int tot = s->num_atoms + s->num_deleted_H;
    printf("N %d %d %d %d %d\n", s->num_atoms, s->num_deleted_H, s->charge, s->bIsotopic, s->ti.num_t_groups);
    for (int i = 0; i < tot; i++) {
        inp_ATOM *a = s->at + i;
        printf("a%d %d %s %d %d %d %d %d %d %d %d %d %d %d %d %d |", i, a->el_number, a->elname, a->orig_at_number, a->orig_compt_at_numb, a->component,
               a->num_H, a->valence, a->chem_bonds_valence, a->at_type, a->endpoint, a->iso_atw_diff, a->num_iso_H[0], a->num_iso_H[1], a->num_iso_H[2], a->p_parity);
        for (int k = 0; k < a->valence; k++) printf(" %d,%d", a->neighbor[k], a->bond_type[k]);
        printf(" |");
        for (int k = 0; k < 4; k++) printf(" %d", a->p_orig_at_num[k]);
        printf(" |");
        for (int k = 0; k < 3; k++) printf(" %d,%d,%d,%d", a->sb_parity[k], a->sb_ord[k], a->sn_ord[k], a->sn_orig_at_num[k]);
        printf("\n");
    }
    if (s->st) for (int i = 0; i < s->num_atoms; i++) {
        inp_ATOM_STEREO *a = s->st + i;
        printf("s%d %d |", i, a->p_parity);
        for (int k = 0; k < 4; k++) printf(" %d", a->p_orig_at_num[k]);
        printf(" |");
        for (int k = 0; k < 3; k++) printf(" %d,%d,%d,%d", a->sb_parity[k], a->sb_ord[k], a->sn_ord[k], a->sn_orig_at_num[k]);
        printf("\n");
    }
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int cs = 0; cs < ncases; cs++) {
        int t;
        printf("#CASE %d\n", cs);
        scanf("%d", &t);
        if (t == 61) {
            int mode; scanf("%d %s", &mode, hx);
            int nn = unhex(hx, buf);
            INCHI_IOSTREAM ios; inchi_ios_init(&ios, INCHI_IOS_TYPE_STRING, NULL);
            ios.s.pStr = buf; ios.s.nUsedLength = nn; ios.s.nPtr = 0;
            SEGM_LINE line; memset(&line, 0, sizeof(line));
            char *pstr = NULL; int state = -77, nerr = -9;
            InpInChI one; memset(&one, 0, sizeof(one));
            int stdf = -5, hso = -5; unsigned char bits = 0;
            int r = InChILine2Data(&ios, &line, &pstr, &state, &nerr, one.pInpInChI, one.nNumComponents, one.nNumProtons, one.s, 1, 1, 0, &stdf, &hso, &bits, &one.polymer, &one.v3000);
            int ok = (r == RI_ERR_EOL || r == RI_ERR_EOF) && !nerr && (one.nNumComponents[0][0] + one.nNumComponents[0][1]);
            printf("R %d\n", ok);
            if (ok) for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) for (int i = 0; i < one.nNumComponents[a][b] && one.pInpInChI[a][b]; i++) {
                INChI *c = &one.pInpInChI[a][b][i];
                if (!c->nNumberOfAtoms) continue;
                INChI *mob = NULL;
                if ((mode & 1) && b == 0 && i < one.nNumComponents[a][1] && one.pInpInChI[a][1] && one.pInpInChI[a][1][i].nNumberOfAtoms) mob = &one.pInpInChI[a][1][i];
                StrFromINChI s; memset(&s, 0, sizeof(s));
                int rr = RestoreAtomConnectionsSetStereo(&s, i, (mode & 2) ? 5 : 0, c, mob);
                printf("G%d%d.%d %d\n", a, b, i, rr);
                if (rr > 0) {
                    int r2 = SetStereoBondTypesFrom0DStereo(&s, c);
                    printf("B %d\n", r2);
                    dump_atoms(&s);
                }
            }
        }
    }
    return 0;
}
"""
HARNESS = r2.HARNESS.split("int main(void)")[0] + MAIN


def main():
    corpus = [json.loads(l) for l in open(ROOT / "tests/data/inchi_corpus.jsonl")]
    inchis = sorted({c["inchi"] for c in corpus if c.get("inchi", "").startswith("InChI=1") and "/z" not in c["inchi"] and "AuxInfo" not in c["inchi"]})
    cases = []
    for ic in inchis:
        for mode in (0, 1, 2, 3):
            cases.append(f"61 {mode} {hx(ic.encode())}")
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        kept, out = [], []
        for c in cases:
            try:
                r = subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True, timeout=5)
            except subprocess.TimeoutExpired:
                continue
            if r.returncode != 0:
                continue
            kept.append(c)
            out += [f"#CASE {len(kept) - 1}"] + r.stdout.decode(errors="replace").split("\n")[1:]
        print(len(cases) - len(kept), "cases dropped")
        cases = kept
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_rvr2.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
