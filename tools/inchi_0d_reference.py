"""Generates tests/data/inchi_0d.json from the InChI C code: the conversion of
0D stereo descriptors into atom/bond parities (`Extract0DParities` of
readinch.c with `FixUnkn0DStereoBonds`, `ReconcileAllCmlBondParities`,
`ReconcileCmlIncidentBondParities` and `get_opposite_sb_atom` of ichister.c)
on random molecules with random, partly invalid descriptors. The function
texts are cut out of the sources. Needs a C compiler and the fetched InChI
sources."""

import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "upstream/external/InChI/INCHI_BASE/src"

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mode.h"
#include "incomdef.h"
#include "ichidrp.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichierr.h"
#include "inchi_api.h"
#include "extr_ct.h"
#ifndef SB_PARITY_FLAG
#define SB_PARITY_FLAG  0x38
#define SB_PARITY_SHFT  3
#define SB_PARITY_MASK  0x07
#endif
#ifndef ATOM_PARITY_WELL_DEF
#define ATOM_PARITY_WELL_DEF(X)     (1 <= (X) && (X) <= 2)
#endif
@FUNCS@
static void hex(const char *s) { printf("x"); for (; s && *s; s++) printf("%02x", (unsigned char) *s); }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, n0d, unk;
        scanf("%d %d %d", &n, &n0d, &unk);
        inp_ATOM *at = calloc(n, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, orig;
            scanf("%d %d %d %d", &v, &cbv, &el, &orig);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].orig_at_number = orig;
            for (int j = 0; j < v; j++) { int x, b; scanf("%d %d", &x, &b); at[i].neighbor[j] = x; at[i].bond_type[j] = b; }
        }
        inchi_Stereo0D *st = calloc(n0d ? n0d : 1, sizeof(inchi_Stereo0D));
        for (int i = 0; i < n0d; i++) {
            int a, b, c2, d, ce, t, p;
            scanf("%d %d %d %d %d %d %d", &a, &b, &c2, &d, &ce, &t, &p);
            st[i].neighbor[0] = a; st[i].neighbor[1] = b; st[i].neighbor[2] = c2; st[i].neighbor[3] = d;
            st[i].central_atom = ce; st[i].type = t; st[i].parity = p;
        }
        char pStrErr[STR_ERR_LEN]; memset(pStrErr, 0, sizeof pStrErr);
        int err = 0;
        Extract0DParities(at, n, st, n0d, pStrErr, &err, unk);
        printf("%d ", err); hex(pStrErr); printf("\n");
        for (int i = 0; i < n; i++) {
            inp_ATOM *a = &at[i];
            printf("A %d %d %d | %d %d %d | %d %d %d | %d %d %d | %d | %d %d %d %d |", a->sb_parity[0], a->sb_parity[1], a->sb_parity[2],
                   a->sb_ord[0], a->sb_ord[1], a->sb_ord[2], a->sn_ord[0], a->sn_ord[1], a->sn_ord[2],
                   a->sn_orig_at_num[0], a->sn_orig_at_num[1], a->sn_orig_at_num[2], a->p_parity,
                   a->p_orig_at_num[0], a->p_orig_at_num[1], a->p_orig_at_num[2], a->p_orig_at_num[3]);
            for (int k = 0; k < a->valence; k++) printf(" %d", a->bond_stereo[k]);
            printf("\n");
        }
        free(at); free(st);
    }
    return 0;
}
"""


def cut(src, start, end_marker="\n}\n"):
    a = src.index(start)
    return src[a:src.index(end_marker, a) + 3]


def build_case(rng):
    n = rng.randrange(3, 10)
    edges = {}
    for i in range(1, n):
        edges[(rng.randrange(i), i)] = 1
    for _ in range(rng.randrange(0, 3)):
        a, b = rng.randrange(n), rng.randrange(n)
        if a != b:
            edges[(min(a, b), max(a, b))] = 1
    for e in edges:
        edges[e] = rng.choice([1, 1, 2, 2, 3, 4])
    adj = [[] for _ in range(n)]
    for (a, b), t in edges.items():
        if len(adj[a]) < 5 and len(adj[b]) < 5:
            adj[a].append((b, t))
            adj[b].append((a, t))
    for l in adj:
        rng.shuffle(l)
    atoms = []
    for i in range(n):
        cbv = sum(t for _, t in adj[i])
        atoms.append({"v": len(adj[i]), "cbv": cbv, "el": 6, "orig": rng.choice([i + 1, i + 1, rng.randrange(1, 40)]), "nb": adj[i]})
    desc = []

    def par():
        return rng.choice([1, 2, 1, 2, 3, 4, 0, 5, 9, 1 | (1 << 3), 2 | (2 << 3), 1 | (4 << 3), 3 | (3 << 3), 2 | (1 << 3), 7, 0x38])

    for (a, b), t in edges.items():
        if t == 2 and (a, b) in [(x, y) for (x, y) in edges]:
            na = [x for x, _ in adj[a] if x != b]
            nb = [x for x, _ in adj[b] if x != a]
            if na and nb and rng.random() < 0.85:
                desc.append([rng.choice(na), a, b, rng.choice(nb), -1, 1, par()])
    for i in range(n):
        if len(adj[i]) >= 3 and rng.random() < 0.5:
            nbs = [x for x, _ in adj[i]]
            rng.shuffle(nbs)
            nbs = nbs[:4]
            while len(nbs) < 4:
                nbs.append(i)
            rng.shuffle(nbs)
            desc.append(nbs + [i, 2, par()])
    # allenes: a=c=b
    for c in range(n):
        dn = [x for x, t in adj[c] if t == 2]
        if len(dn) == 2 and rng.random() < 0.9:
            a, b = dn
            na = [x for x, _ in adj[a] if x != c]
            nb = [x for x, _ in adj[b] if x != c]
            if na and nb:
                desc.append([rng.choice(na), a, b, rng.choice(nb), c, 3, par()])
    # corruption
    for d in desc:
        if rng.random() < 0.12:
            k = rng.randrange(4)
            d[k] = rng.choice([-1, n, 0, n - 1, 100])
        if rng.random() < 0.05:
            d[4] = rng.choice([-1, 0, n])
        if rng.random() < 0.05:
            d[5] = rng.choice([0, 1, 2, 3, 9])
    rng.shuffle(desc)
    return {"n": n, "atoms": atoms, "desc": desc, "unk": rng.choice([3, 4])}


def main():
    rng = random.Random(41)
    cases = [build_case(rng) for _ in range(1500)]
    # tetrahedral descriptors with a central atom outside the structure are undefined behaviour in C
    for c in cases:
        c["desc"] = [d for d in c["desc"] if not (d[5] == 2 and not (0 <= d[4] < c["n"]))]
    stdin = f"{len(cases)}\n"
    for c in cases:
        stdin += f"{c['n']} {len(c['desc'])} {c['unk']}\n"
        for a in c["atoms"]:
            stdin += f"{a['v']} {a['cbv']} {a['el']} {a['orig']} " + " ".join(f"{x} {t}" for x, t in a["nb"]) + "\n"
        for d in c["desc"]:
            stdin += " ".join(map(str, d)) + "\n"
    ri = (SRC / "readinch.c").read_text(errors="replace")
    ic = (SRC / "ichister.c").read_text(errors="replace")
    start = ri.index("int Extract0DParities( inp_ATOM *at,", ri.index("inchi_Stereo0D * CreateInchi_Stereo0D"))
    funcs = [ri[start:ri.index("\n}\n", start) + 3]]
    funcs.append(cut(ic, "int FixUnkn0DStereoBonds( inp_ATOM *at, int num_at )\n{"))
    funcs.append(cut(ic, "int ReconcileAllCmlBondParities( inp_ATOM *at,\n                                 int num_atoms,\n                                 int bDisconnected )\n{"))
    funcs.append(cut(ic, "int ReconcileCmlIncidentBondParities( inp_ATOM *at,\n                                      int cur_atom,\n                                      int prev_atom,\n                                      S_CHAR *visited,\n                                      int bDisconnected )\n{"))
    funcs.append(cut(ic, "int get_opposite_sb_atom( inp_ATOM *at,\n                          int cur_atom,\n                          int icur2nxt,\n                          int *pnxt_atom,\n                          int *pinxt2cur,\n                          int *pinxt_sb_parity_ord )\n{"))
    protos = """int FixUnkn0DStereoBonds( inp_ATOM *at, int num_at );
int ReconcileAllCmlBondParities( inp_ATOM *at, int num_atoms, int bDisconnected );
int ReconcileCmlIncidentBondParities( inp_ATOM *at, int cur_atom, int prev_atom, S_CHAR *visited, int bDisconnected );
int get_opposite_sb_atom( inp_ATOM *at, int cur_atom, int icur2nxt, int *pnxt_atom, int *pinxt2cur, int *pinxt_sb_parity_ord );
int Extract0DParities( inp_ATOM *at, int nNumAtoms, inchi_Stereo0D *stereo0D, int num_stereo0D, char *pStrErr, int *err, int vABParityUnknown );
"""
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS.replace("@FUNCS@", protos + "\n".join(funcs)))
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    k = 0
    data = []
    for c in cases:
        head = out[k]
        k += 1
        atoms = out[k:k + c["n"]]
        k += c["n"]
        data.append({"in": c, "head": head, "atoms": atoms})
    path = ROOT / "tests/data/inchi_0d.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
