"""Generates tests/data/inchi_onestruct.json from the InChI C code:
`CreateOneStructureINChI` of runichi.c (preprocessing, component extraction,
`CreateOneComponentINChI` per component, reuse of the disconnected-layer results
for the reconnected layer) on the random salts, metal complexes, ion pairs and
molecules of `inchi_preprocess.json`. Links the upstream InChI base library
(see inchi_markalt_reference.py).

    python tools/inchi_onestruct_reference.py
"""
import importlib.util
import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


cr = load("create", "tools/inchi_create_reference.py")
ma = cr.ma
SRC = ma.SRC

HARNESS = cr.HARNESS.split("int main(void)")[0] + r"""
#include "ichimain.h"
#include "ichierr.h"
#include "ichi_io.h"
int CreateOneStructureINChI( CANON_GLOBALS *pCG, INCHI_CLOCK *ic, STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle,
                             PINChI2 *pINChI2[INCHI_NUM], PINChI_Aux2 *pINChI_Aux2[INCHI_NUM], int iINChI,
                             INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file,
                             ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, COMP_ATOM_DATA composite_norm_data2[][TAUT_NUM + 1],
                             long num_inp, INCHI_IOS_STRING *strbuf, NORM_CANON_FLAGS *pncFlags );
int main(void) {
    int ncases;
    setvbuf(stdout, NULL, _IONBF, 0);
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, nbonds;
        unsigned long long tf; int fnu, nowarn;
        unsigned long long mode; int cmpc;
        scanf("%d %d %llu %d %d %llu %d", &n, &nbonds, &tf, &fnu, &nowarn, &mode, &cmpc);
        inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, chg, rad, iso, nh, ih0, ih1, ih2;
            double x, y, z; char en[8];
            scanf("%d %d %d %d %d %d %d %d %d %d %lf %lf %lf %s", &v, &cbv, &el, &chg, &rad, &iso, &nh, &ih0, &ih1, &ih2, &x, &y, &z, en);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].charge = chg;
            at[i].radical = rad; at[i].iso_atw_diff = iso; at[i].num_H = nh;
            at[i].num_iso_H[0] = ih0; at[i].num_iso_H[1] = ih1; at[i].num_iso_H[2] = ih2;
            at[i].x = x; at[i].y = y; at[i].z = z; at[i].orig_at_number = i + 1;
            strcpy(at[i].elname, en);
            for (int j = 0; j < v; j++) { int q, b, s; scanf("%d %d %d", &q, &b, &s); at[i].neighbor[j] = q; at[i].bond_type[j] = b; at[i].bond_stereo[j] = s; }
        }
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        cg.m_num_bit = 16;
        cg.m_bBit = calloc(16, sizeof(bitWord));
        for (int i = 0; i < 16; i++) cg.m_bBit[i] = 1 << i;
        INCHI_CLOCK ic; memset(&ic, 0, sizeof(ic));
        INPUT_PARMS ip; memset(&ip, 0, sizeof(ip));
        ip.nMode = mode; ip.bTautFlags = tf; ip.bFixNonUniformDraw = fnu; ip.bNoWarnings = nowarn; ip.bCompareComponents = cmpc;
        STRUCT_DATA sd; memset(&sd, 0, sizeof(sd));
        ORIG_ATOM_DATA orig; memset(&orig, 0, sizeof(orig));
        ORIG_ATOM_DATA prep[2]; memset(prep, 0, sizeof(prep));
        orig.at = at; orig.num_inp_atoms = n; orig.num_inp_bonds = nbonds;
        PINChI2 *pi[INCHI_NUM]; PINChI_Aux2 *pa[INCHI_NUM];
        memset(pi, 0, sizeof(pi)); memset(pa, 0, sizeof(pa));
        COMP_ATOM_DATA comp[INCHI_NUM][TAUT_NUM + 1]; memset(comp, 0, sizeof(comp));
        NORM_CANON_FLAGS nc; memset(&nc, 0, sizeof(nc));
        char title[256]; title[0] = 0;
        INCHI_IOSTREAM s1, s2, s3, s4;
        inchi_ios_init(&s1, INCHI_IOS_TYPE_STRING, NULL); inchi_ios_init(&s2, INCHI_IOS_TYPE_STRING, NULL);
        inchi_ios_init(&s3, INCHI_IOS_TYPE_STRING, NULL); inchi_ios_init(&s4, INCHI_IOS_TYPE_STRING, NULL);
        printf("#CASE %d\n", c);
        for (int ii = 0; ii < 2; ii++) {
            int ret = CreateOneStructureINChI(&cg, &ic, &sd, &ip, title, pi, pa, ii, &s1, &s2, &s3, &s4, &orig, prep, comp, 1, NULL, &nc);
            printf("S%d %d %d %d %d [%s] %d %d %d %d %d\n", ii, ret, sd.nErrorCode, sd.nErrorType, sd.nStructReadError, sd.pStrErrStruct,
                   sd.num_non_taut[ii], sd.num_taut[ii], sd.num_components[ii], prep[0].num_components, prep[1].num_components);
            for (int j = 0; j < sd.num_components[ii] && pi[ii]; j++) {
                for (int k = 0; k < 2; k++) {
                    INChI *p = pi[ii][j][k]; INChI_Aux *a = pa[ii][j][k];
                    if (!p) { printf("I%d.%d.%d -\n", ii, j, k); continue; }
                    printf("I%d.%d.%d %d %d %d %s %d %d %d %d", ii, j, k, p->nErrorCode, p->nNumberOfAtoms, p->nTotalCharge, p->szHillFormula ? p->szHillFormula : "NULL",
                           p->lenConnTable, p->lenTautomer, p->nNumberOfIsotopicAtoms, p->nNumberOfIsotopicTGroups);
                    printf(" |"); for (int q = 0; q < p->lenConnTable; q++) printf(" %d", p->nConnTable[q]);
                    printf(" |"); for (int q = 0; q < p->lenTautomer; q++) printf(" %d", p->nTautomer[q]);
                    printf(" |"); for (int q = 0; q < p->nNumberOfAtoms; q++) printf(" %d", p->nNum_H[q]);
                    printf(" | %d %d\n", a ? a->nNumberOfAtoms : -1, a ? a->nNumRemovedProtons : -1);
                }
            }
        }
        printf("F"); for (int ii = 0; ii < 2; ii++) for (int k = 0; k < 2; k++) printf(" %llu/%llu/%llu", (unsigned long long)nc.bTautFlags[ii][k], (unsigned long long)nc.bTautFlagsDone[ii][k], (unsigned long long)nc.bNormalizationFlags[ii][k]);
        printf("\n");
    }
    return 0;
}
"""


def fmt(c, mode, cmpc):
    atoms = c["atoms"]
    nbonds = c["nbonds"]
    s = f"{len(atoms)} {nbonds} {c['flags']} {c['fnu']} {c['nowarn']} {mode} {cmpc}\n"
    for a in atoms:
        s += (f"{len(a['nb'])} {a['cbv']} {a['el']} {a['chg']} {a['rad']} {a['iso']} {a['nh']} "
              f"{a['ih'][0]} {a['ih'][1]} {a['ih'][2]} {a['xyz'][0]} {a['xyz'][1]} {a['xyz'][2]} {a['sym']} "
              + " ".join(f"{x} {t} {st}" for x, t, st in a["nb"]) + "\n")
    return s


def main():
    rng = random.Random(31337)
    base = json.load(open(ROOT / "tests/data/inchi_preprocess.json"))
    cases = []
    for rec in base:
        c = rec["in"]
        if not c["atoms"]:
            continue
        cases.append((c, rng.choice(cr.MODES), rng.choice([0, 0, 1])))
    objs = ma.build_lib()
    inputs = [fmt(*c) for c in cases]
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{ma.SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        keep, lines = [], []
        r = subprocess.run([str(exe)], input=(f"{len(inputs)}\n" + "".join(inputs)).encode(), capture_output=True)
        if r.returncode == 0:
            out = r.stdout.decode().strip().split("\n")
            keep = list(range(len(inputs)))
        else:
            out = []
            for i, t in enumerate(inputs):
                r1 = subprocess.run([str(exe)], input=("1\n" + t).encode(), capture_output=True)
                if r1.returncode == 0:
                    keep.append(i)
                    out += [f"#CASE {len(keep) - 1}"] + r1.stdout.decode().strip().split("\n")[1:]
            print(len(inputs) - len(keep), "cases dropped")
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(keep)
    data = [{"in": inputs[i], "out": out[a + 1:b]} for i, a, b in zip(keep, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_onestruct.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
