#pragma once
// Vertex 1 (production): gamma(beam) + p1(LH2 target, at rest)
//                           -> p(bachelor) + Lambda + Lambdabar

#include "core/fourvector.h"
#include "core/rng.h"

namespace mcgen::rescatter {

struct Vertex1Result {
    FourVec proton;
    FourVec lambda;
    FourVec lambdabar;
};

// beamEnergy: photon energy (GeV) along +z. Flat 3-body phase space via
// core::kin3body; returns false if beamEnergy is below the pLambdaLambdabar
// threshold (should not happen for beam energies in the documented 8-9 GeV
// working range).
bool generateVertex1(Rng& rng, double beamEnergy, Vertex1Result& out);

}  // namespace mcgen::rescatter
