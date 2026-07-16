// Physics sanity tests for the two-vertex rescattering generator:
//  - global 4-momentum conservation for gamma + p1(rest) + p2(rest) ->
//    p Lambda K+ pi+ pi- pi0 (the single most important check per
//    context/02-generator-design.md in the pantilamb workspace: this is
//    what makes MM^2(gamma+p1->detected) peak at m_p^2, not 0)
//  - each final-state particle is on its mass shell
// Exits nonzero on any failure.
#include <cmath>
#include <cstdio>

#include "core/rng.h"
#include "core/targetgeom.h"
#include "rescatter/event.h"
#include "rescatter/particles.h"

namespace {
int failures = 0;

void checkClose(double a, double b, double tol, const char* what) {
    if (std::abs(a - b) > tol) {
        std::printf("FAIL: %s (%.9g vs %.9g)\n", what, a, b);
        ++failures;
    }
}

void checkOnShell(const mcgen::FourVec& p, double mass, const char* what) {
    checkClose(p.e * p.e - p.p2(), mass * mass, 1e-6, what);
}
}  // namespace

int main() {
    using namespace mcgen;
    using namespace mcgen::rescatter;

    Rng rng(20260716);
    TargetGeom geom;
    geom.icode = 2;
    geom.size1 = 2.0;
    geom.size2 = 2.0;
    geom.size3 = 51.0;
    geom.size4 = 79.0;
    geom.jcode = 1;
    geom.xsigma = 0.2;
    geom.ysigma = 0.2;

    int generated = 0;
    for (int i = 0; i < 5000; ++i) {
        Event ev;
        Truth truth;
        if (!generateRescatterEvent(rng, geom, 8.0, 9.0, ev, truth)) continue;
        ++generated;

        // Global initial state: gamma(+z) + two protons at rest.
        FourVec lhs{0.0, 0.0, truth.beamEnergy, 0.0,
                    truth.beamEnergy + 2.0 * kMassProton};

        FourVec rhs{};
        for (size_t k = 1; k < ev.size(); ++k) {
            rhs.px += ev[k].p.px;
            rhs.py += ev[k].p.py;
            rhs.pz += ev[k].p.pz;
            rhs.e += ev[k].p.e;
        }
        checkClose(rhs.px, lhs.px, 1e-6, "global px");
        checkClose(rhs.py, lhs.py, 1e-6, "global py");
        checkClose(rhs.pz, lhs.pz, 1e-6, "global pz");
        checkClose(rhs.e, lhs.e, 1e-6, "global E");

        checkOnShell(ev[1].p, kMassProton, "proton on-shell");
        checkOnShell(ev[2].p, kMassLambda, "Lambda on-shell");
        checkOnShell(ev[3].p, kMassKaonCh, "K+ on-shell");
        checkOnShell(ev[4].p, kMassPionCh, "pi+ on-shell");
        checkOnShell(ev[5].p, kMassPionCh, "pi- on-shell");
        checkOnShell(ev[6].p, kMassPion0, "pi0 on-shell");
        checkOnShell(truth.lambdabar, kMassLambda, "Lambdabar on-shell (truth)");
    }

    if (generated < 4900)
        std::printf("WARN: only %d/5000 attempts succeeded\n", generated);
    if (failures == 0) std::printf("all rescatter tests passed (%d events)\n",
                                   generated);
    return failures == 0 ? 0 : 1;
}
