#pragma once
// Vertex 2 (annihilation): Lambdabar + p2(second LH2 target proton, at rest)
//                             -> K+ pi+ pi- pi0
//
// Kinematically isomorphic to vertex 1's own beam+target -> N-body problem:
// a moving projectile (here Lambdabar, produced at vertex 1) hits a proton
// at rest. Reuses core::kin45body exactly as-is.

#include "core/fourvector.h"
#include "core/rng.h"

namespace mcgen::rescatter {

struct Vertex2Result {
    FourVec kaonPlus;
    FourVec piPlus;
    FourVec piMinus;
    FourVec pi0;
};

// lambdabar: the Lambdabar four-vector produced at vertex 1 (lab frame).
// Flat 4-body phase space (GENBOD) boosted into the lab frame of
// lambdabar + p2. Returns false if below the K+pi+pi-pi0 threshold (should
// not happen: available energy here is set by E_Lambdabar, ~2.3-3 GeV,
// comfortably above threshold).
bool generateVertex2(Rng& rng, const FourVec& lambdabar, Vertex2Result& out);

}  // namespace mcgen::rescatter
