// Driver for tests/data/coordgen_minim.txt: CoordgenMinimizer::minimizeMolecule on a molecule with perturbed coordinates
// (after forceUpdateStruct and CoordgenFragmenter::splitIntoFragments); prints the interaction counts, energies and the
// final coordinates.
//   c++ -std=c++17 -O1 -ffp-contract=off -I SRC minim.cpp OBJDIR/*.o -o minim
#define private public
#define protected public
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include "CoordgenFragmenter.h"
#include "CoordgenMinimizer.h"
#include "sketcherMinimizerFragment.h"
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
            int a, b, o;
            sscanf(line.c_str(), "B %d %d %d", &a, &b, &o);
            auto* bd = mol.addNewBond(mol._atoms[a], mol._atoms[b]);
            bd->setBondOrder(o);
        }
        std::getline(std::cin, line); // END
        sketcherMinimizerMolecule::forceUpdateStruct(mol._atoms, mol._bonds, mol._rings);
        CoordgenFragmenter::splitIntoFragments(&mol);
        CoordgenMinimizer min;
        min.m_evenAngles = false;
        min.m_atoms = mol._atoms;
        min.m_bonds = mol._bonds;
        min.m_molecules.push_back(&mol);
        min.minimizeMolecule(&mol);
        float e = min.scoreInteractions();
        printf("interactions %zu clash %zu stretch %zu bend %zu energy %08x\n", min._interactions.size(),
               min._intramolecularClashInteractions.size(), min._stretchInteractions.size(), min._bendInteractions.size(), bits(e));
        for (size_t i = 0; i < mol._atoms.size(); i++)
            printf("a %zu %08x %08x\n", i, bits(mol._atoms[i]->coordinates.x()), bits(mol._atoms[i]->coordinates.y()));
        printf("$$\n");
    }
    return 0;
}
