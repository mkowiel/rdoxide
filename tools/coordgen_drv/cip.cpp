// Driver for tests/data/coordgen_cip.txt: CIP priorities, stereo reading/writing and double bond first neighbours.
//   c++ -std=c++17 -O1 -ffp-contract=off -I SRC cip.cpp OBJDIR/*.o -o cip
// Input per molecule: "MOL n m", "A Z charge xbits ybits", "B a b order wedge" (wedge 0 none, 1 wedge, 2 dashed, from the
// first atom), "END".
#define private public
#define protected public
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerBond.h"
#include "sketcherMinimizerMolecule.h"
#include "sketcherMinimizerRing.h"

static float fromBits(uint32_t u)
{
    float f;
    memcpy(&f, &u, 4);
    return f;
}

int main()
{
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.rfind("MOL", 0) != 0) continue;
        int n, m;
        sscanf(line.c_str(), "MOL %d %d", &n, &m);
        sketcherMinimizerMolecule mol;
        for (int i = 0; i < n; i++) {
            std::getline(std::cin, line);
            int z, ch;
            unsigned xb, yb;
            sscanf(line.c_str(), "A %d %d %x %x", &z, &ch, &xb, &yb);
            auto* a = mol.addNewAtom();
            a->setAtomicNumber(z);
            a->charge = ch;
            a->setCoordinates(sketcherMinimizerPointF(fromBits(xb), fromBits(yb)));
        }
        for (int i = 0; i < m; i++) {
            std::getline(std::cin, line);
            int a, b, o, w;
            sscanf(line.c_str(), "B %d %d %d %d", &a, &b, &o, &w);
            auto* bd = mol.addNewBond(mol._atoms[a], mol._atoms[b]);
            bd->setBondOrder(o);
            if (w) {
                bd->hasStereochemistryDisplay = true;
                bd->isWedge = w == 1;
            }
        }
        std::getline(std::cin, line); // END
        sketcherMinimizerMolecule::forceUpdateStruct(mol._atoms, mol._bonds, mol._rings);
        auto& A = mol._atoms;
        auto& B = mol._bonds;
        auto ai = [&](sketcherMinimizerAtom* a) -> long { return a ? std::find(A.begin(), A.end(), a) - A.begin() : -1; };
        for (size_t i = 0; i < A.size(); i++) {
            auto& nb = A[i]->neighbors;
            for (size_t j = 0; j + 1 < nb.size(); j++)
                for (size_t k = j + 1; k < nb.size(); k++)
                    printf("p %zu %zu %zu %ld\n", i, j, k, ai(sketcherMinimizerAtom::CIPPriority(nb[j], nb[k], A[i])));
        }
        for (size_t i = 0; i < A.size(); i++) {
            auto* a = A[i];
            int r = a->readStereochemistry(false);
            printf("r %zu %d stereogenic %d ring %d inv %d rs", i, r, a->m_isStereogenic, a->_hasRingChirality, a->m_clockwiseInvert);
            for (int p : a->m_RSPriorities) printf(" %d", p);
            printf("\n");
        }
        for (size_t i = 0; i < A.size(); i++) {
            auto* a = A[i];
            a->hasStereochemistrySet = true;
            a->isR = i % 2 == 0;
            a->writeStereoChemistry();
            printf("w %zu set %d ring %d\n", i, a->hasStereochemistrySet, a->_hasRingChirality);
        }
        for (size_t i = 0; i < B.size(); i++) {
            auto* b = B[i];
            printf("b %zu disp %d wedge %d rev %d first %ld %ld\n", i, b->hasStereochemistryDisplay, b->isWedge, b->isReversed,
                   ai(b->startAtomCIPFirstNeighbor()), ai(b->endAtomCIPFirstNeighbor()));
        }
        printf("$$\n");
    }
    return 0;
}
