"""Generates tests/data/inchi_rvr1.json from the InChI C code: the first part
of the InChI -> structure restoration (ichirvr1.c): `GetAtomRestoreInfo`
(with `clean_charge_val`, `get_sp_element_type`) on random small atom
graphs, and `GetTgroupInfoFromInChI` / `FillOutpStructEndpointFromInChI`
on the InChI strings of the official RDKit corpus. The harness includes
ichiread.c. Links the upstream InChI base library (see
inchi_markalt_reference.py).

    python tools/inchi_rvr1_reference.py
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


r3 = load("read3", "tools/inchi_read3_reference.py")
r2 = r3.r2
ma = r3.ma
SRC = ma.SRC

PRE = r3.MAIN.split("int main(void)")[0]
MAIN = PRE + r"""
static void dump_va(VAL_AT *v) {
    printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", v->cDoNotAddH, v->cMetal, v->cNumBondsToMetal, v->cInitFlowToMetal,
           v->cInitValenceToMetal, v->cInitOrigValenceToMetal, v->cMaxFlowToMetal, v->cInitFreeValences, v->cInitCharge, v->cNumValenceElectrons,
           v->cPeriodicRowNumber, v->cMinRingSize, v->cPeriodicNumber, v->cnListIndex, v->nCMinusGroupEdge, v->nCPlusGroupEdge, v->nMetalGroupEdge, v->nTautGroupEdge);
}
static void dump_ti(T_GROUP_INFO *ti) {
    printf("TI %d %d %d | ", ti->num_t_groups, ti->nNumEndpoints, ti->max_num_t_groups);
    for (int i = 0; i < ti->num_t_groups && ti->t_group; i++)
        printf("(%d %d %d %d %d) ", ti->t_group[i].num[0], ti->t_group[i].num[1], ti->t_group[i].nGroupNumber, ti->t_group[i].nNumEndpoints, ti->t_group[i].nFirstEndpointAtNoPos);
    printf("| ");
    for (int i = 0; i < ti->nNumEndpoints && ti->nEndpointAtomNumber; i++) printf("%d,", ti->nEndpointAtomNumber[i]);
    printf("| ");
    if (ti->tGroupNumber) for (int i = 0; i < ti->num_t_groups; i++) printf("%d,%d ", ti->tGroupNumber[i], ti->tGroupNumber[TGSO_SYMM_IORDER * ti->num_t_groups + i]); else printf("NULL");
    printf("\n");
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int cs = 0; cs < ncases; cs++) {
        int t;
        printf("#CASE %d\n", cs);
        scanf("%d", &t);
        if (t == 20) {
            int mob, n; scanf("%d %d", &mob, &n);
            inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
            AT_NUMB *ep = calloc(n + 1, sizeof(AT_NUMB));
            for (int i = 0; i < n; i++) {
                int el, nh, cbv, e, nb;
                scanf("%d %d %d %d %d", &el, &nh, &cbv, &e, &nb);
                at[i].el_number = el; at[i].num_H = nh; at[i].chem_bonds_valence = cbv; at[i].endpoint = e; ep[i] = e; at[i].valence = nb;
                for (int k = 0; k < nb; k++) { int nn, bt; scanf("%d %d", &nn, &bt); at[i].neighbor[k] = nn; at[i].bond_type[k] = bt; }
            }
            VAL_AT *va = calloc(n + 1, sizeof(VAL_AT));
            SRM srm; SetUpSrm(&srm);
            CANON_GLOBALS cg; memset(&cg, 0, sizeof(cg));
            for (int i = 0; i < n; i++) {
                int j;
                va[i].cNumValenceElectrons = get_sp_element_type(at[i].el_number, &j);
                va[i].cPeriodicRowNumber = j;
                va[i].cPeriodicNumber = at[i].el_number;
                va[i].cNumValenceElectrons--;
                if (is_el_a_metal(at[i].el_number)) va[i].cMetal = 1;
            }
            for (int i = 0; i < n; i++) {
                va[i].cInitFreeValences = 0;
                int r = GetAtomRestoreInfo(&cg, at, i, va, &srm, mob, ep);
                printf("A%d %d %d | ", i, r, at[i].chem_bonds_valence); dump_va(&va[i]); printf("\n");
            }
        } else if (t == 21) {
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
                int na = c->nNumberOfAtoms;
                inp_ATOM *at = calloc(na + 2, sizeof(inp_ATOM));
                AT_NUMB *ep = calloc(na + 2, sizeof(AT_NUMB));
                T_GROUP_INFO ti; memset(&ti, 0, sizeof(ti));
                for (int pass = 0; pass < 2; pass++) {
                    int rr = GetTgroupInfoFromInChI(&ti, mode & 1 ? at : NULL, mode & 2 ? ep : NULL, c);
                    printf("G%d%d.%d.%d %d ", a, b, i, pass, rr); dump_ti(&ti);
                    printf("E "); for (int k = 0; k < na; k++) printf("%d,", at[k].endpoint); printf(" | ");
                    for (int k = 0; k < na; k++) printf("%d,", ep[k]); printf("\n");
                }
                AT_NUMB *ep2 = NULL;
                int rr = FillOutpStructEndpointFromInChI(c, &ep2);
                printf("F %d ", rr); for (int k = 0; k < na; k++) printf("%d,", ep2[k]); printf("\n");
                free(ep2);
            }
        } else if (t == 22) {
            for (int e = 0; e < 130; e++) { int row; int ty = get_sp_element_type(e, &row); printf("%d,%d ", ty, row); }
            printf("\n");
        }
    }
    return 0;
}
"""
HARNESS = r2.HARNESS.split("int main(void)")[0] + MAIN

R = random.Random(8642)
hx = r2.hx

ELS = [1, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 14, 15, 15, 16, 16, 17, 34, 35, 53, 11, 12, 26, 29, 13, 3, 2, 33, 52, 20, 30, 50, 82, 79, 92]


def atom_case():
    n = R.randrange(1, 7)
    mob = R.choice([0, 1])
    adj = [[] for _ in range(n)]
    # random tree + extra bonds
    edges = {}
    for i in range(1, n):
        j = R.randrange(0, i)
        edges[(j, i)] = R.choice([1, 1, 1, 2, 2, 3, 4])
    for _ in range(R.randrange(0, 3)):
        i, j = sorted(R.sample(range(n), 2)) if n > 1 else (0, 0)
        if i != j:
            edges[(i, j)] = R.choice([1, 1, 2, 4])
    for (i, j), b in edges.items():
        adj[i].append((j, b))
        adj[j].append((i, b))
    parts = [f"20 {mob} {n}"]
    for i in range(n):
        el = R.choice(ELS)
        nh = R.choice([0, 0, 1, 2, 3])
        s = sum(min(b, 3) if b != 4 else 1 for _, b in adj[i])
        cbv = s + R.choice([0, 0, 0, 1, 2]) if any(b == 4 for _, b in adj[i]) else s + R.choice([0, 0, 0, 1])
        ep = R.choice([0, 0, 0, 1, 2])
        parts.append(f"{el} {nh} {cbv} {ep} {len(adj[i])} " + " ".join(f"{j} {b}" for j, b in adj[i]))
    return " ".join(parts)


def main():
    corpus = [json.loads(l) for l in open(ROOT / "tests/data/inchi_corpus.jsonl")]
    inchis = sorted({c["inchi"] for c in corpus if c.get("inchi", "").startswith("InChI=1") and "/z" not in c["inchi"] and "AuxInfo" not in c["inchi"]})
    cases = [atom_case() for _ in range(40000)]
    for ic in inchis:
        for mode in (3, 0, 1, 2):
            cases.append(f"21 {mode} {hx(ic.encode())}")
    cases.append("22")
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        kept, out = [], []
        for c in cases:
            try:
                r1 = subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True, timeout=2)
            except subprocess.TimeoutExpired:
                continue
            if r1.returncode != 0:
                continue
            kept.append(c)
            out += [f"#CASE {len(kept) - 1}"] + r1.stdout.decode(errors="replace").split("\n")[1:]
        print(len(cases) - len(kept), "cases dropped")
        cases = kept
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_rvr1.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
