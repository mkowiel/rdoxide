// Driver for tests/data/pubchem_shape_align.txt: pubchem-align3d (NCBI, daefab3) on two sets of atoms.
//   c++ -std=c++17 -O1 -ffp-contract=off -I SRC align.cpp SRC/shape_functions1.cpp SRC/shape_functions2.cpp SRC/shape_neighbor.cpp -o align
// Input per case: "CASE nref nfit color opt_hex pre post", then nref + nfit lines "xbits ybits zbits radius_hex type"
// (coordinates as float bit patterns). Output: the shape/colour volumes of both sets, the qrat of the steric frames,
// and the result of Neighbor_Conformers (matrix as float bits, nbr_st, nbr_ct as double bits).
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "shape_functions.hpp"

static double dfromhex(const char *s) { uint64_t u = strtoull(s, nullptr, 16); double d; memcpy(&d, &u, 8); return d; }
static uint64_t dbits(double d) { uint64_t u; memcpy(&u, &d, 8); return u; }
static uint32_t fbits(float f) { uint32_t u; memcpy(&u, &f, 4); return u; }

struct Side {
  std::vector<float> coord; std::vector<double> rad, alpha; std::vector<unsigned> types, vol;
  std::map<unsigned, std::vector<unsigned>> cmap; double sov = 0, sof = 0;
};

static void read(Side &s, int n) {
  char line[256];
  for (int i = 0; i < n; i++) {
    std::cin.getline(line, sizeof line);
    unsigned xb, yb, zb, t; char rh[40];
    sscanf(line, "%x %x %x %39s %u", &xb, &yb, &zb, rh, &t);
    float x, y, z; memcpy(&x, &xb, 4); memcpy(&y, &yb, 4); memcpy(&z, &zb, 4);
    s.coord.push_back(x); s.coord.push_back(y); s.coord.push_back(z);
    s.rad.push_back(dfromhex(rh)); s.types.push_back(t);
  }
  Align3D::setAlpha(s.rad.data(), n, s.alpha);
  Align3D::getVolumeAtomIndexVector(s.types.data(), n, s.vol);
  Align3D::getColorAtomType2IndexVectorMap(s.types.data(), n, s.cmap);
  s.sov = Align3D::ComputeShapeOverlap(s.coord.data(), s.alpha, s.vol, s.coord.data(), s.alpha, s.vol);
  s.sof = Align3D::ComputeFeatureOverlap(s.coord.data(), s.alpha, s.cmap, s.coord.data(), s.alpha, s.cmap);
}

int main() {
  Align3D::setUseCutOff(true);
  std::string line;
  while (std::getline(std::cin, line)) {
    if (line.rfind("CASE", 0) != 0) continue;
    int nref, nfit, color, pre, post; char oh[40];
    sscanf(line.c_str(), "CASE %d %d %d %39s %d %d", &nref, &nfit, &color, oh, &pre, &post);
    double opt = dfromhex(oh);
    Side r, f; read(r, nref); read(f, nfit);
    std::set<unsigned> joint; Align3D::getJointColorTypeSet(r.types.data(), nref, f.types.data(), nfit, joint);
    auto rm = r.cmap, fm = f.cmap;
    Align3D::restrictColorAtomType2IndexVectorMap(rm, joint);
    Align3D::restrictColorAtomType2IndexVectorMap(fm, joint);
    unsigned rq, fq;
    { auto c = r.coord; Align3D::TransformCoordToStericFrameAndCalculateQrat(nref, 0.7225f, c.data(), rq); }
    { auto c = f.coord; Align3D::TransformCoordToStericFrameAndCalculateQrat(nfit, 0.7225f, c.data(), fq); }
    float matrix[12]; double st, ct;
    Align3D::Neighbor_Conformers(r.coord.data(), r.alpha, r.vol, rm, r.sov, r.sof, f.coord.data(), f.alpha, f.vol, fm,
                                 f.sov, f.sof, color && !joint.empty(), pre, post, opt, matrix, st, ct);
    printf("vol %016llx %016llx %016llx %016llx\n", (unsigned long long)dbits(r.sov), (unsigned long long)dbits(r.sof),
           (unsigned long long)dbits(f.sov), (unsigned long long)dbits(f.sof));
    printf("qrat %u %u\n", rq, fq);
    printf("m");
    for (int i = 0; i < 12; i++) printf(" %08x", fbits(matrix[i]));
    printf("\nst %016llx ct %016llx\n$$\n", (unsigned long long)dbits(st), (unsigned long long)dbits(ct));
  }
}
