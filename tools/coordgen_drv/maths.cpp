// Driver for tests/data/coordgen_maths.txt: calls sketcherMinimizerMaths.h of coordgenlibs 3.0.2 on pseudo random input.
//   c++ -std=c++17 -O1 -ffp-contract=off -I coordgenlibs-3.0.2 tools/coordgen_drv/maths.cpp -o maths
// Every line: <name> <inputs as hex float bits> ; <outputs as hex float bits>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include "sketcherMinimizerMaths.h"

static uint64_t state = 12345;
static float rnd(float scale = 100.f)
{
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    int v = static_cast<int>((state >> 33) % 20001) - 10000;
    return v * scale / 10000.f;
}
static uint32_t bits(float f)
{
    uint32_t u;
    memcpy(&u, &f, 4);
    return u;
}
static void out(const char* name, const std::vector<float>& in, const std::vector<float>& res, int b = -1)
{
    printf("%s", name);
    for (float f : in) printf(" %08x", bits(f));
    printf(" ;");
    for (float f : res) printf(" %08x", bits(f));
    if (b >= 0) printf(" %d", b);
    printf("\n");
}
using P = sketcherMinimizerPointF;
using M = sketcherMinimizerMaths;
static P rp() { float x = rnd(), y = rnd(); return P(x, y); }

int main()
{
    for (int i = 0; i < 300; i++) {
        float f = rnd(300.f);
        out("round2", {f}, {roundToTwoDecimalDigits(f)});
        for (int p : {0, 1, 3, 4}) {
            std::vector<float> in = {f, static_cast<float>(p)};
            out("roundp", in, {roundToPrecision(f, p)});
        }
        P a = rp(), b = rp(), c = rp(), d = rp();
        {
            P t = a;
            t.normalize();
            out("normalize", {a.x(), a.y()}, {a.length(), t.x(), t.y()});
            float s = rnd(1.f), co = rnd(1.f);
            t = a;
            t.rotate(s, co);
            out("rotate", {a.x(), a.y(), s, co}, {t.x(), t.y()});
            out("parallel", {a.x(), a.y(), b.x(), b.y()}, {P(a).parallelComponent(b).x(), P(a).parallelComponent(b).y()});
            t = a;
            t.round();
            P t3 = a;
            t3.round(3);
            out("pround", {a.x(), a.y()}, {t.x(), t.y(), t3.x(), t3.y()});
        }
        {
            P r;
            bool hit = M::intersectionOfSegments(a, b, c, d, &r);
            out("intersect", {a.x(), a.y(), b.x(), b.y(), c.x(), c.y(), d.x(), d.y()}, {r.x(), r.y()}, hit);
        }
        out("signedAngle", {a.x(), a.y(), b.x(), b.y(), c.x(), c.y()}, {M::signedAngle(a, b, c)});
        out("unsignedAngle", {a.x(), a.y(), b.x(), b.y(), c.x(), c.y()}, {M::unsignedAngle(a, b, c)});
        out("coincide", {a.x(), a.y(), a.x() + 0.00001f, a.y()}, {}, M::pointsCoincide(a, P(a.x() + 0.00001f, a.y())));
        out("sameSide", {a.x(), a.y(), b.x(), b.y(), c.x(), c.y(), d.x(), d.y()}, {}, M::sameSide(a, b, c, d));
        {
            P pr = M::projectPointOnLine(a, b, c);
            out("project", {a.x(), a.y(), b.x(), b.y(), c.x(), c.y()}, {pr.x(), pr.y()});
            float t;
            float sd = M::squaredDistancePointSegment(a, b, c, &t);
            out("sqdistseg", {a.x(), a.y(), b.x(), b.y(), c.x(), c.y()}, {sd, t});
            out("sqdist", {a.x(), a.y(), b.x(), b.y()}, {M::squaredDistance(a, b)});
            P mp = M::mirrorPoint(a, b, c);
            out("mirror", {a.x(), a.y(), b.x(), b.y(), c.x(), c.y()}, {mp.x(), mp.y()});
            float tt = (rnd(1.f) + 1.f) / 2.f;
            P bz = M::pointOnCubicBezier(a, b, c, d, tt);
            out("bezier", {a.x(), a.y(), b.x(), b.y(), c.x(), c.y(), d.x(), d.y(), tt}, {bz.x(), bz.y()});
        }
        out("dot", {}, {}, 0);
        {
            float x1 = rnd(), y1 = rnd(), z1 = rnd(), x2 = rnd(), y2 = rnd(), z2 = rnd(), x3 = rnd(), y3 = rnd(), z3 = rnd(),
                  x4 = rnd(), y4 = rnd(), z4 = rnd();
            out("length3D", {x1, y1, z1}, {M::length3D(x1, y1, z1)});
            out("dot3D", {x1, y1, z1, x2, y2, z2}, {M::dotProduct3D(x1, y1, z1, x2, y2, z2)});
            float xr, yr, zr;
            M::crossProduct3D(x1, y1, z1, x2, y2, z2, xr, yr, zr);
            out("cross3D", {x1, y1, z1, x2, y2, z2}, {xr, yr, zr});
            out("dist3D", {x1, y1, z1, x2, y2, z2}, {M::distance3D(x1, y1, z1, x2, y2, z2)});
            out("angle3D", {x1, y1, z1, x2, y2, z2, x3, y3, z3}, {M::angle3D(x1, y1, z1, x2, y2, z2, x3, y3, z3)});
            out("dihedral3D", {x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4},
                {M::dihedral3D(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)});
            float dl = rnd(), dm = rnd(), dn = rnd();
            float dd = M::length3D(dl, dm, dn);
            if (dd > 0.f) { dl /= dd; dm /= dd; dn /= dd; }
            out("cannon", {x1 / 10, y1 / 10, z1 / 10, dl, dm, dn, x2 / 10, y2 / 10, z2 / 10, 1.f, 1.5f},
                {M::cannonBallDistance(x1 / 10, y1 / 10, z1 / 10, dl, dm, dn, x2 / 10, y2 / 10, z2 / 10, 1.f, 1.5f)});
        }
    }
    for (int n = 3; n < 14; n++) {
        for (int rep = 0; rep < 6; rep++) {
            std::vector<P> knots;
            std::vector<float> in;
            for (int i = 0; i < n; i++) { knots.push_back(rp()); in.push_back(knots.back().x()); in.push_back(knots.back().y()); }
            std::vector<P> f, s;
            M::ClosedBezierControlPoints(knots, f, s);
            std::vector<float> res;
            for (auto& p : f) { res.push_back(p.x()); res.push_back(p.y()); }
            for (auto& p : s) { res.push_back(p.x()); res.push_back(p.y()); }
            out("closedbezier", in, res);
            std::vector<float> a(n), b(n), c(n), rhs(n);
            for (int i = 0; i < n; i++) { a[i] = rnd(1.f); b[i] = 4.f + rnd(1.f); c[i] = rnd(1.f); rhs[i] = rnd(); }
            std::vector<float> in2;
            for (int i = 0; i < n; i++) { in2.push_back(a[i]); in2.push_back(b[i]); in2.push_back(c[i]); in2.push_back(rhs[i]); }
            out("tridiag", in2, M::tridiagonalSolve(a, b, c, rhs));
            out("cyclic", in2, M::cyclicSolve(a, b, c, 1.f, 1.f, rhs));
        }
    }
    return 0;
}
