"""Generates tests/data/inchi_strings.json from the InChI C code: the string
helpers of util.c (`normalize_string`, `lrtrim`, `read_upto_delim`,
`extract_inchi_substring`, `inchi_stricmp`, ..., `MakeRemovedProtonsString`)
on random strings. Strings are passed as hex. Needs a C compiler and the
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
#include <string.h>
#include "mode.h"
#include "util.h"
#include "ichierr.h"
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
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        char h1[512], h2[512], s1[512] = {0}, s2[512] = {0}, b[1024];
        int slen, maxlen, nrp, bi;
        int e[3], r[3];
        scanf("%s %s %d %d %d %d %d %d %d %d %d %d", h1, h2, &slen, &maxlen, &nrp, &e[0], &e[1], &e[2], &r[0], &r[1], &r[2], &bi);
        memset(s1, 0, sizeof s1); memset(s2, 0, sizeof s2); unhex(h1, s1); unhex(h2, s2);
        int len;
        strcpy(b, s1); len = normalize_string(b); hex(b, strlen(b)); printf(" %d|", len);
        strcpy(b, s1); len = dotify_non_printable_chars(b); hex(b, strlen(b)); printf(" %d|", len);
        strcpy(b, s1); lrtrim(b, &len); hex(b, strlen(b)); printf(" %d|", len);
        strcpy(b, s1); remove_trailing_spaces(b); hex(b, strlen(b)); printf("|");
        strcpy(b, s1); remove_one_lf(b); hex(b, strlen(b)); printf("|");
        {
            char *p = s1; char field[512] = "";
            int n = read_upto_delim(&p, field, maxlen, ", \t");
            printf("%d ", n); hex(field, n > 0 ? n : 0); printf(" %d|", p ? (int)(p - s1) : -1);
        }
        {
            char *o = NULL;
            extract_inchi_substring(&o, s1, slen);
            if (o) hex(o, strlen(o)); else printf("N");
            printf("|");
            o = NULL;
            extract_auxinfo_substring(&o, s1, slen);
            if (o) hex(o, strlen(o)); else printf("N");
            printf("|");
        }
        {
            int a = inchi_stricmp(s1, s2), m = inchi_memicmp(s1, s2, strlen(s1) < strlen(s2) ? strlen(s1) : strlen(s2));
            printf("%d %d|", a < 0 ? -1 : a > 0, m < 0 ? -1 : m > 0);
        }
        {
            char t[64]; memset(t, 'z', 63); t[63] = 0;
            mystrncpy(t, s1, maxlen);
            hex(t, 64); printf("|");
        }
        {
            char out[512]; int niso = 0;
            NUM_H ex[3] = {e[0], e[1], e[2]}, rm[3] = {r[0], r[1], r[2]};
            int l = MakeRemovedProtonsString(nrp, (bi & 2) ? NULL : ex, (bi & 4) ? NULL : rm, bi & 1, out, &niso);
            hex(out, l); printf(" %d", niso);
        }
        printf("\n");
    }
    return 0;
}
"""

ALPHABET = list(" \t\t\n\r\x0b\x0cabcAZ09()=/;,%@?InChI=AuxInfo=") + ["InChI=1S/CH4/h1H4", "AuxInfo=1/0/N:1", " ", "  ", "\x01", "\x7f", "\xe9"]


def rstr(rng):
    out = ""
    for _ in range(rng.randrange(0, 14)):
        out += rng.choice(ALPHABET)
    if rng.random() < 0.3:
        out = rng.choice(["", " ", "\n"]) + "InChI=1S/C2H6/c1-2/h1-2H3" + rng.choice(["", " x", "\r\n", "%", " \t"]) + out
    if rng.random() < 0.15:
        out += "AuxInfo=1/1/N:1,2 rest"
    return out


def hx(s):
    return "x" + s.encode("latin-1").hex()


def main():
    rng = random.Random(3)
    cases = []
    for _ in range(600):
        s1, s2 = rstr(rng), rstr(rng)
        if rng.random() < 0.3:
            s2 = s1.swapcase()
        cases.append({
            "s1": s1, "s2": s2, "slen": rng.choice([0, 3, 10, 40, 200]), "maxlen": rng.choice([1, 3, 10, 64]),
            "nrp": rng.choice([0, 0, 2, -3]), "e": [rng.choice([0, 0, 1, 2]) for _ in range(3)],
            "r": [rng.choice([0, 0, 1, 3]) for _ in range(3)], "bi": rng.randrange(8)})
    stdin = f"{len(cases)}\n"
    for c in cases:
        stdin += f"{hx(c['s1'])} {hx(c['s2'])} {c['slen']} {c['maxlen']} {c['nrp']} " + " ".join(map(str, c["e"] + c["r"])) + f" {c['bi']}\n"
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        # MakeRemovedProtonsString is compiled out under COMPILE_ANSI_ONLY
        # (set for TARGET_API_LIB), so take its text from util.c
        util = (SRC / "util.c").read_text(errors="replace")
        a = util.index("int MakeRemovedProtonsString(")
        fn = util[a:util.index("#endif", a)]
        (Path(d) / "h.c").write_text(HARNESS.replace("int main(void) {", fn + "\nint main(void) {", 1))
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode("ascii"), capture_output=True, check=True).stdout.decode("latin-1").strip().split("\n")
    for c, line in zip(cases, out):
        c["out"] = line.split("|")
        c["s1"] = hx(c["s1"])
        c["s2"] = hx(c["s2"])
    path = ROOT / "tests/data/inchi_strings.json"
    path.write_text(json.dumps(cases, separators=(",", ":")))
    print(len(cases), "cases ->", path)


main()
