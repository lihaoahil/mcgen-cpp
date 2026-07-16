#include "rescatter/event.h"

#include "rescatter/particles.h"
#include "rescatter/vertex1.h"
#include "rescatter/vertex2.h"

namespace mcgen::rescatter {
namespace {

// 1-based indexing, matching Generator's own convention (generator.cpp:539,
// 548-549): the pseudo-particle occupies slot 1, so its direct children
// carry iparent=1, not 0. GEN2HDDM_new indexes parent-1 into its particle
// table, so an iparent of 0 underflows -- confirmed by an out_of_range
// crash when this was first tried with iparent=0.
Particle makeFinal(const FourVec& p, int lund, const Vec3& r, int index) {
    Particle q;
    q.p = p;
    q.lund = lund;
    q.rCreate = r;
    q.pDecay = p;   // no generator-level in-flight decay applied (v1)
    q.rDecay = r;
    q.index = index;
    q.iparent = 1;  // direct child of the pseudo-particle at slot 1
    return q;
}

}  // namespace

bool generateRescatterEvent(Rng& rng, const TargetGeom& geom, double beamLo,
                            double beamHi, Event& outEvent, Truth& outTruth) {
    const double beamEnergy = rng.uniform(beamLo, beamHi);

    Vertex1Result v1;
    if (!generateVertex1(rng, beamEnergy, v1)) return false;

    // v1 simplification: vertex 2 co-located with vertex 1 (see STATUS.md).
    Vertex2Result v2;
    if (!generateVertex2(rng, v1.lambdabar, v2)) return false;

    const Vec3 r0 = getVertex(rng, geom);

    outEvent.clear();
    outEvent.push_back(Particle{});
    outEvent[0].lund = 999;
    outEvent[0].index = 1;
    outEvent[0].rCreate = r0;
    outEvent[0].rDecay = r0;
    outEvent.push_back(makeFinal(v1.proton, kLundProton, r0, 2));
    outEvent.push_back(makeFinal(v1.lambda, kLundLambda, r0, 3));
    outEvent.push_back(makeFinal(v2.kaonPlus, kLundKaonPlus, r0, 4));
    outEvent.push_back(makeFinal(v2.piPlus, kLundPiPlus, r0, 5));
    outEvent.push_back(makeFinal(v2.piMinus, kLundPiMinus, r0, 6));
    outEvent.push_back(makeFinal(v2.pi0, kLundPi0, r0, 7));

    outTruth.lambdabar = v1.lambdabar;
    outTruth.vertex1 = r0;
    outTruth.vertex2 = r0;
    outTruth.beamEnergy = beamEnergy;
    return true;
}

}  // namespace mcgen::rescatter
