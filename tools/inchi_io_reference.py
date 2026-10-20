"""Generates tests/data/inchi_io.json from the InChI C code: the string
stream readers (`inchi_ios_gets*`, `inchi_sgets`, `inchi_fgetsLf`,
`inchi_strbuf_addline`), the print buffers and `Output_RecordInfo` on random
inputs. Texts are passed as hex. Needs a C compiler and the fetched InChI
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
#include "ichi_io.h"
#include "util.h"
int inchi_ios_str_getc(INCHI_IOSTREAM *ios);
const char gsMissing[] = "is missing";
const char gsEmpty[] = "";
const char gsSpace[] = " ";
const char gsEqual[] = "=";
static void hex(const char *s, int n) {
    printf("x");
    for (int i = 0; i < n; i++) printf("%02x", (unsigned char) s[i]);
}
static int unhex(const char *h, char *out) {
    int n = 0;
    for (h++; h[0] && h[1]; h += 2) { unsigned v; sscanf(h, "%2x", &v); out[n++] = (char) v; }
    out[n] = 0;
    return n;
}
static void show_ios(INCHI_IOSTREAM *s) {
    printf(" [%d %d %d]", s->s.nAllocatedLength, s->s.nUsedLength, s->s.nPtr);
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        char h[2048], text[2048];
        char kind[8];
        scanf("%s", kind);
        if (kind[0] == 'R') {
            /* record info */
            int num, nolbl; unsigned long id;
            char hl[256], hv[256], lf[256], tab[256], l[256], v[256];
            int hasl, hasv;
            scanf("%d %d %lu %d %s %d %s", &num, &nolbl, &id, &hasl, hl, &hasv, hv);
            unhex(hl, l); unhex(hv, v);
            INCHI_IOSTREAM o; inchi_ios_init(&o, INCHI_IOS_TYPE_STRING, NULL);
            Output_RecordInfo(&o, num, nolbl, hasl ? l : NULL, hasv ? v : NULL, id, "\n", "\t");
            hex(o.s.pStr ? o.s.pStr : "", o.s.pStr ? strlen(o.s.pStr) : 0); show_ios(&o);
            printf("\n");
            continue;
        }
        scanf("%s", h);
        int n = unhex(h, text);
        if (kind[0] == 'P') {
            /* print sequence: pieces separated by length-prefixed parts; text split in chunks of 7 */
            INCHI_IOSTREAM o; inchi_ios_init(&o, INCHI_IOS_TYPE_STRING, NULL);
            INCHI_IOS_STRING sb; inchi_strbuf_init(&sb, 16, 8);
            char chunk[16];
            for (int i = 0; i < n; i += 7) {
                int l = n - i < 7 ? n - i : 7;
                memcpy(chunk, text + i, l); chunk[l] = 0;
                inchi_ios_print(&o, "%s", chunk);
                inchi_strbuf_printf(&sb, "%s", chunk);
                if (i == 14) inchi_strbuf_printf_from(&sb, 3, "%s", "ab");
            }
            printf("%d ", o.s.nAllocatedLength); hex(o.s.pStr ? o.s.pStr : "", o.s.pStr ? strlen(o.s.pStr) : 0); show_ios(&o);
            printf(" | %d %d %d ", sb.nAllocatedLength, sb.nUsedLength, sb.nPtr); hex(sb.pStr, strlen(sb.pStr));
            inchi_strbuf_update(&sb, 40); printf(" %d", sb.nAllocatedLength);
            inchi_strbuf_reset(&sb); printf(" %d %d %d", sb.nAllocatedLength, sb.nUsedLength, sb.nPtr);
            inchi_ios_flush_not_displayed(&o);
            printf(" | "); hex(o.s.pStr ? o.s.pStr : "", o.s.pStr ? strlen(o.s.pStr) : 0); show_ios(&o);
            printf("\n");
            continue;
        }
        /* reader sequences */
        INCHI_IOSTREAM in; inchi_ios_init(&in, INCHI_IOS_TYPE_STRING, NULL);
        inchi_ios_print(&in, "%s", text);
        in.s.nPtr = 0;
        INCHI_IOS_STRING sb; inchi_strbuf_init(&sb, 8, 4);
        int nops; scanf("%d", &nops);
        for (int k = 0; k < nops; k++) {
            char op[4]; int a, b2;
            scanf("%s %d %d", op, &a, &b2);
            char line[1024]; int tl = 0; int r;
            switch (op[0]) {
            case 'g': r = inchi_ios_gets(line, a, &in, &tl); printf("%d %d ", r, tl); hex(line, r < 0 ? 0 : strlen(line)); break;
            case 't': r = inchi_ios_getsTab(line, a, &in, &tl); printf("%d %d ", r, tl); hex(line, r < 0 ? 0 : strlen(line)); break;
            case '1': r = inchi_ios_getsTab1(line, a, &in, &tl); printf("%d %d ", r, tl); hex(line, r < 0 ? 0 : strlen(line)); break;
            case 's': { char *p = inchi_sgets(line, a, &in); if (p) { printf("S "); hex(p, strlen(p)); } else printf("N"); break; }
            case 'f': { char *p = inchi_fgetsLf(line, a, &in); if (p) { printf("S "); hex(line, strlen(line)); } else printf("N"); break; }
            case 'c': r = inchi_ios_str_getc(&in); printf("%d", r); break;
            case 'a': r = inchi_strbuf_addline(&sb, &in, a, b2); printf("%d ", r); hex(sb.pStr, strlen(sb.pStr)); printf(" %d", sb.nAllocatedLength); break;
            }
            printf(" #%d;", in.s.nPtr);
        }
        printf("\n");
    }
    return 0;
}
"""


def hx(s):
    return "x" + s.encode("latin-1").hex()


def text(rng):
    pieces = ["abc", " ", "\t", "\n", "\r", "\r\n", "x y", "  ", "InChI=1S/C", "0123456789", "z"]
    return "".join(rng.choice(pieces) for _ in range(rng.randrange(0, 18)))


def main():
    rng = random.Random(21)
    cases = []
    lines = []
    for _ in range(400):
        t = text(rng)
        ops = []
        for _ in range(rng.randrange(1, 8)):
            op = rng.choice("gt1sfcaa")
            if op in "gt1":
                ops.append((op, rng.choice([4, 6, 10, 40]), 0))
            elif op in "sf":
                ops.append((op, rng.choice([1, 2, 5, 12, 100]), 0))
            elif op == "a":
                ops.append((op, rng.randrange(2), rng.randrange(2)))
            else:
                ops.append((op, 0, 0))
        cases.append(("Q", t, ops))
        lines.append(f"Q {hx(t)} {len(ops)} " + " ".join(f"{o} {a} {b}" for o, a, b in ops))
    for _ in range(100):
        t = text(rng) or "q"
        cases.append(("P", t, None))
        lines.append(f"P {hx(t)}")
    for _ in range(100):
        lab = rng.choice(["", "ID", "name"])
        val = rng.choice(["", "v1", "7"])
        num = rng.randrange(1, 50)
        nol = rng.choice([0, 0, 0, 1])
        sid = rng.choice([0, 0, 12, 345])
        cases.append(("R", [lab, val, num, nol, sid], None))
        lines.append(f"R {num} {nol} {sid} {1 if lab else 0} {hx(lab)} {1 if val else 0} {hx(val)}")
    stdin = f"{len(lines)}\n" + "\n".join(lines) + "\n"
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS)
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "ichi_io.c"), str(SRC / "util.c"), str(SRC / "ichierr.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode("ascii"), capture_output=True, check=True).stdout.decode("latin-1").strip().split("\n")
    data = []
    for (k, t, ops), o in zip(cases, out):
        if k == "Q":
            data.append({"kind": k, "text": hx(t), "ops": [list(x) for x in ops], "out": o})
        elif k == "P":
            data.append({"kind": k, "text": hx(t), "out": o})
        else:
            data.append({"kind": k, "args": t, "out": o})
    path = ROOT / "tests/data/inchi_io.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
