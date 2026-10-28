// Driver for tests/data/coordgen_sketch.txt: what RDKit's CoordGen::addCoords does with coordgenlibs (External/CoordGen/
// CoordGen.h): a molecule from atoms, bonds and double bond stereo, sketcherMinimizer(0.01) initialize and
// runGenerateCoordinates; prints the clean pose flag and the coordinates.
//   c++ -std=c++17 -O1 -ffp-contract=off -I SRC sketch.cpp OBJDIR/*.o -o sketch
// Input per molecule: "MOL n m", "A Z charge xbits ybits", "B a b order", "S bond atom1 atom2 c|t", "END".
#define private public
#define protected public
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include "sketcherMinimizer.h"
#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerBond.h"
#include "sketcherMinimizerMolecule.h"

static uint32_t bits(float f)
{
    uint32_t u;
    memcpy(&u, &f, 4);
    return u;
}

int main()
{
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.rfind("MOL", 0) != 0) continue;
        int n, m;
        sscanf(line.c_str(), "MOL %d %d", &n, &m);
        sketcherMinimizer minimizer(0.01f);
        auto* min_mol = new sketcherMinimizerMolecule();
        minimizer.setTreatNonterminalBondsToMetalAsZOBs(false);
        std::vector<sketcherMinimizerAtom*> ats(n);
        std::vector<sketcherMinimizerBond*> bnds(m);
        for (int i = 0; i < n; i++) {
            std::getline(std::cin, line);
            int z, ch;
            unsigned xb, yb;
            sscanf(line.c_str(), "A %d %d %x %x", &z, &ch, &xb, &yb);
            auto* a = min_mol->addNewAtom();
            a->molecule = min_mol;
            a->atomicNumber = z;
            a->charge = ch;
            ats[i] = a;
        }
        for (int i = 0; i < m; i++) {
            std::getline(std::cin, line);
            int a, b, o;
            sscanf(line.c_str(), "B %d %d %d", &a, &b, &o);
            auto* bd = min_mol->addNewBond(ats[a], ats[b]);
            bd->bondOrder = o;
            bnds[i] = bd;
        }
        min_mol->assignBondsAndNeighbors(ats, bnds);
        std::getline(std::cin, line);
        while (line.rfind("S ", 0) == 0) {
            int bi, a1, a2;
            char c;
            sscanf(line.c_str(), "S %d %d %d %c", &bi, &a1, &a2, &c);
            sketcherMinimizerBondStereoInfo sinfo;
            sinfo.atom1 = ats[a1];
            sinfo.atom2 = ats[a2];
            sinfo.stereo = c == 'c' ? sketcherMinimizerBondStereoInfo::cis : sketcherMinimizerBondStereoInfo::trans;
            bnds[bi]->setStereoChemistry(sinfo);
            bnds[bi]->setAbsoluteStereoFromStereoInfo();
            std::getline(std::cin, line);
        }
        minimizer.initialize(min_mol);
        bool clean = minimizer.runGenerateCoordinates();
        printf("clean %d\n", clean);
        for (int i = 0; i < n; i++) printf("a %d %08x %08x\n", i, bits(ats[i]->coordinates.x()), bits(ats[i]->coordinates.y()));
        printf("$$\n");
    }
    return 0;
}
