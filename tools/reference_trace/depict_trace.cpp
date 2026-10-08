// Re-runs RDDepict::compute2DCoords (no coordMap/templates/sampling) with the
// official RDKit libraries, printing exact coordinates after each stage.
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/EmbeddedFrag.h>
#include <GraphMol/Depictor/DepictUtils.h>
#include <GraphMol/Rings.h>
#include <GraphMol/MolOps.h>
#include <cstdio>
#include <list>
using namespace RDDepict;

static void dump(const char *stage, std::list<EmbeddedFrag> &efrags) {
  int fi = 0;
  for (auto &f : efrags) {
    for (auto &ea : f.GetEmbeddedAtoms()) {
      printf("%s f%d %u %.17g %.17g\n", stage, fi, ea.first, ea.second.loc.x, ea.second.loc.y);
    }
    ++fi;
  }
}

int main(int argc, char **argv) {
  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(argv[1]));
  RDKit::ROMol cp(*mol);
  std::vector<int> atomRanks(cp.getNumAtoms());
  for (auto i = 0u; i < cp.getNumAtoms(); ++i) atomRanks[i] = getAtomDepictRank(cp.getAtomWithIdx(i));
  RDKit::VECT_INT_VECT arings;
  RDKit::MolOps::symmetrizeSSSR(cp, arings, true);
  RDKit::MolOps::assignStereochemistry(cp, false);
  std::list<EmbeddedFrag> efrags;
  if (!arings.empty()) {
    RDKit::INT_INT_VECT_MAP neighMap;
    RingUtils::makeRingNeighborMap(arings, neighMap);
    auto cnrs = arings.size();
    boost::dynamic_bitset<> fusDone(cnrs);
    auto curr = 0u;
    while (curr < cnrs) {
      RDKit::INT_VECT fused;
      RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
      RDKit::VECT_INT_VECT frings;
      for (auto rid : fused) frings.push_back(arings.at(rid));
      EmbeddedFrag efrag(&cp, frings, false);
      efrag.setupNewNeighs();
      efrags.push_back(efrag);
      size_t rix;
      for (rix = 0; rix < cnrs; ++rix) if (!fusDone[rix]) { curr = rix; break; }
      if (rix == cnrs) break;
    }
  }
  for (auto bond : cp.bonds()) {
    if (bond->getBondType() == RDKit::Bond::DOUBLE && bond->getStereo() > RDKit::Bond::STEREOANY &&
        !cp.getRingInfo()->numBondRings(bond->getIdx()) && bond->getStereoAtoms().size() == 2) {
      EmbeddedFrag efrag(bond);
      efrag.setupNewNeighs();
      efrags.push_back(efrag);
    }
  }
  RDKit::INT_LIST nratms;
  boost::dynamic_bitset<> done(cp.getNumAtoms());
  for (auto &f : efrags) for (auto &ea : f.GetEmbeddedAtoms()) done[ea.first] = 1;
  for (auto aid = 0u; aid < cp.getNumAtoms(); ++aid) if (!done[aid]) nratms.push_back(aid);
  auto largest = [&]() {
    std::list<EmbeddedFrag>::iterator mfri = efrags.end(); int msiz = 0;
    for (auto e = efrags.begin(); e != efrags.end(); ++e) if (!e->isDone() && e->Size() > msiz) { msiz = e->Size(); mfri = e; }
    return mfri;
  };
  auto mri = largest();
  while (mri != efrags.end() || nratms.size() > 0) {
    if (mri == efrags.end()) {
      int mrank = RDKit::MAX_INT; auto mnri = nratms.end();
      for (auto nri = nratms.begin(); nri != nratms.end(); ++nri) {
        int rank = atomRanks.at(*nri) * cp.getNumAtoms() + *nri;
        if (rank < mrank) { mrank = rank; mnri = nri; }
      }
      EmbeddedFrag efrag((*mnri), &cp);
      nratms.erase(mnri);
      efrags.push_back(efrag);
      mri = efrags.end(); --mri;
    }
    mri->markDone();
    mri->expandEfrag(nratms, efrags);
    mri = largest();
  }
  dump("initial", efrags);
  for (auto &e : efrags) e.removeCollisionsBondAndSpiroFlip();
  dump("bondflip", efrags);
  for (auto &e : efrags) {
    e.removeCollisionsOpenAngles();
    dump("openangles", efrags);
    e.removeCollisionsShortenBonds();
    dump("shorten", efrags);
    e.removeCollisionsPathAngleExpansion();
    dump("pathangle", efrags);
  }
  return 0;
}
