"""Generates tests/data/inchi_ring.json from the InChI C code: the smallest
ring through each bond (`is_bond_in_Nmax_memb_ring`) and the 3-membered ring
test (`is_atom_in_3memb_ring`) for random graphs. Needs a C compiler and the
fetched InChI sources."""

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
#include "mode.h"
#include "ichiring.h"
#include "inpdef.h"
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n;
        scanf("%d", &n);
        inp_ATOM *at = calloc(n, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, rs, na;
            scanf("%d %d %d", &v, &rs, &na);
            at[i].valence = v;
            at[i].nRingSystem = rs;
            at[i].nNumAtInRingSystem = na;
            for (int j = 0; j < v; j++) { int x; scanf("%d", &x); at[i].neighbor[j] = x; }
        }
        QUEUE *q = QueueCreate(n + 1, sizeof(AT_RANK));
        AT_RANK *lev = calloc(n, sizeof(AT_RANK));
        S_CHAR *src = calloc(n, 1);
        int nq;
        scanf("%d", &nq);
        for (int k = 0; k < nq; k++) {
            int a, o, m;
            scanf("%d %d %d", &a, &o, &m);
            printf("%d ", is_bond_in_Nmax_memb_ring(at, a, o, q, lev, src, m));
        }
        printf("|");
        for (int i = 0; i < n; i++) printf(" %d", is_atom_in_3memb_ring(at, i));
        printf("\n");
        QueueDelete(q);
        free(at); free(lev); free(src);
    }
    return 0;
}
"""


def graph(rng):
    n = rng.randrange(2, 14)
    adj = [[] for _ in range(n)]
    # random tree plus extra edges
    for i in range(1, n):
        j = rng.randrange(i)
        adj[i].append(j)
        adj[j].append(i)
    for _ in range(rng.randrange(0, n)):
        a, b = rng.randrange(n), rng.randrange(n)
        if a != b and b not in adj[a] and len(adj[a]) < 4 and len(adj[b]) < 4:
            adj[a].append(b)
            adj[b].append(a)
    for a in adj:
        rng.shuffle(a)
    return adj


def main():
    rng = random.Random(11)
    cases = []
    for _ in range(400):
        adj = graph(rng)
        n = len(adj)
        rs = [rng.randrange(1, 3) for _ in range(n)]
        nar = [rng.randrange(0, 6) for _ in range(n)]
        qs = []
        for a in range(n):
            for o in range(len(adj[a])):
                qs.append([a, o, rng.choice([0, 2, 3, 4, 5, 6, 8, 20])])
        cases.append({"adj": adj, "rs": rs, "nar": nar, "q": qs})
    stdin = f"{len(cases)}\n"
    for c in cases:
        stdin += f"{len(c['adj'])}\n"
        for i, a in enumerate(c["adj"]):
            stdin += f"{len(a)} {c['rs'][i]} {c['nar'][i]} " + " ".join(map(str, a)) + "\n"
        stdin += f"{len(c['q'])}\n" + "\n".join(" ".join(map(str, q)) for q in c["q"]) + "\n"
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS)
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "ichiring.c"), str(SRC / "util.c"), str(SRC / "ichierr.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin, capture_output=True, text=True, check=True).stdout.strip().split("\n")
    for c, line in zip(cases, out):
        r, t = line.split("|")
        c["ring"] = [int(x) for x in r.split()]
        c["tri"] = [int(x) for x in t.split()]
    path = ROOT / "tests/data/inchi_ring.json"
    path.write_text(json.dumps(cases, separators=(",", ":")))
    print(len(cases), "graphs ->", path)


main()
