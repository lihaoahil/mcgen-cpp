#pragma once
// Assembles the two vertices into an observable mcgen::Event plus an MC
// truth sidecar carrying what the ASCII stream can't: the (consumed)
// Lambdabar four-vector and the true production/annihilation vertices.

#include "core/fourvector.h"
#include "core/generator.h"  // mcgen::Event / mcgen::Particle
#include "core/rng.h"
#include "core/targetgeom.h"

namespace mcgen::rescatter {

struct Truth {
    FourVec lambdabar;  // true Lambdabar 4-vector (never itself final-state)
    Vec3 vertex1{};      // production vertex (cm)
    Vec3 vertex2{};      // annihilation vertex (cm) -- == vertex1 in v1
                         // (co-located; displaced-vertex sampling deferred,
                         // see STATUS.md known v1 simplifications)
    double beamEnergy = 0.0;
};

// Runs vertex 1 then vertex 2 and assembles the result:
//   outEvent[0]    = pseudo-particle placeholder (lund=999, not written out)
//   outEvent[1..6] = p, Lambda, K+, pi+, pi-, pi0 (the observable final state)
// Returns false if either vertex's phase-space call fails (should not
// happen for beam energies in the documented working range) -- caller
// should just retry.
bool generateRescatterEvent(Rng& rng, const TargetGeom& geom, double beamLo,
                            double beamHi, Event& outEvent, Truth& outTruth);

}  // namespace mcgen::rescatter
