#pragma once
// Breit-Wigner lineshape samplers ported from bwran.F and
// relativisticbwran.F.
//
// bwran: the Fortran sampled 1/(x^2+1) on [-30,30] by inverting a 1001-bin
// cumulative table with bin blur; here the truncated Cauchy is inverted
// analytically (statistically equivalent, no binning artifacts).
//
// RelativisticBWSampler: table inversion is kept (the mass-dependent-width
// shape has no closed-form inverse).  The Fortran carried two copy-pasted
// subroutines (relativisticbwran1/2) because cached tables lived in static
// common blocks; a class instance per resonance replaces that.

#include <string>
#include <vector>

#include "core/rng.h"

namespace mcgen {

// Draw x in [-30,30] from 1/(x^2+1).  Convert to a mass via
// m = (Gamma/2)*x + m_pole, as the Fortran callers do.
double bwran(Rng& rng);

class RelativisticBWSampler {
  public:
    // rmnot: central mass; gammanot: nominal full width Gamma (GeV);
    // rm1, rm2: decay product masses; lvalue: orbital L of the decay.
    RelativisticBWSampler(double rmnot, double gammanot, double rm1,
                          double rm2, int lvalue);

    // Draw a mass from the lineshape (RELATIVISTICBWRAN's bwmass).
    double sample(Rng& rng) const;

    bool belowThreshold() const { return ignoreMasses_; }

  private:
    static constexpr int kNbin = 1000;
    double rm1_, rm2_;
    double rmmin_ = 0.0, drm_ = 0.0;
    bool ignoreMasses_ = false;  // central mass below decay threshold
    std::vector<double> cdf_;    // normalized partial sums, cdf_[0] = 0
};

// RELATIVISTICLINESHAPE: relativistic BW times the two-body phase-space
// factors of the production (resonance + recoil rm3 at total energy ww) and
// of the decay.  Because ww differs event by event the shape cannot be
// pre-tabulated; rejection sampling with an adaptive envelope is used, as in
// the Fortran.  A width >= 9998 requests the pure phase-space shape.
class RelativisticLineshape {
  public:
    struct Result {
        double mass = 0.0;
        int ifail = 0;  // 0 ok; 1 empty range; 2/3/4 kinematics failures
    };

    // rmnot/gammanot/rm1/rm2/lvalue as in RelativisticBWSampler; rm3 is the
    // recoil mass and ww the event's total available energy.
    Result sample(Rng& rng, double rmnot, double gammanot, double rm1,
                  double rm2, int lvalue, double rm3, double ww);

  private:
    double bwmaxval_ = 1000.0;  // adaptive rejection envelope
    double scale_ = 1000.0;
    bool first_ = true;
};

}  // namespace mcgen
