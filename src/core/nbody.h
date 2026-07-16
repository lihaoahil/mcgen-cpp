#pragma once
// N-body phase-space decays, ported from kin3body.F, kin3bodyqf.F,
// kin3bodyhulthen.F, and kin45body.F (which embeds CERNLIB GENBOD).

#include <array>
#include <vector>

#include "core/fourvector.h"
#include "core/rng.h"

namespace mcgen {

// KIN3BODY: flat three-body phase space via Dalitz-variable rejection.
// p0 is the decaying system in the lab; m1..m3 the product masses.
bool kin3body(Rng& rng, const FourVec& p0, double m1, double m2, double m3,
              FourVec& p1, FourVec& p2, FourVec& p3);

// Result bookkeeping the Fortran passed back through the /fermi/ common.
struct QuasiFreeResult {
    double prec = 0.0;    // recoil (spectator) momentum magnitude
    double qfmass = 0.0;  // off-shell mass of the decaying quasi-system
};

// KIN3BODYQF: quasi-free two-body decay with the third particle recoiling
// with Fermi momentum drawn uniformly inside a sphere of radius pfermi.
bool kin3bodyqf(Rng& rng, const FourVec& p0, double m1, double m2, double m3,
                double pfermi, FourVec& p1, FourVec& p2, FourVec& p3,
                QuasiFreeResult& qf);

// KIN3BODYHULTHEN: same, but the recoil momentum follows the deuteron
// Hulthen momentum distribution (sampled below 0.40 GeV/c).
bool kin3bodyhulthen(Rng& rng, const FourVec& p0, double m1, double m2,
                     double m3, FourVec& p1, FourVec& p2, FourVec& p3,
                     QuasiFreeResult& qf);

// GENBOD (Raubold-Lynch, constant-cross-section mode KGENEV=1): decay of a
// system of invariant mass tecm at rest into n bodies with the given masses.
// Fills pcm with (px,py,pz,m,E) vectors; returns the event weight, or a
// negative value when the masses exceed tecm.
double genbod(Rng& rng, double tecm, const std::vector<double>& masses,
              std::vector<FourVec>& pcm);

// KIN45BODY: 4- or 5-body phase space of the moving system p0, boosting the
// GENBOD rest-frame result to the lab.  out.size() selects the multiplicity.
bool kin45body(Rng& rng, const FourVec& p0, const std::vector<double>& masses,
               std::vector<FourVec>& out);

}  // namespace mcgen
