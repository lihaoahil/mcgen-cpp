#pragma once
// Production/decay dynamics, ported from mc_gen.F subroutines: REGGEPICK,
// DOUBLEREGGEPICK, EXPONENTIALPICK, SCULPTMASS, PHASESPACEWEIGHT(3),
// WEAKDECAY, POLARIZE, PRODUCT_STUFFIT, PRODUCT_ROTATE, BOOST_TO_LAB,
// GETT, INVARMASS, GSMEAR.

#include "core/fourvector.h"
#include "core/rng.h"

namespace mcgen {

// REGGEPICK: choose cos(theta_cm) of the produced object so that the t
// distribution follows exp(2 b ln(s/s0) t') above |tcutoff|, with a linear
// rise below it.  Returns false on kinematics failure (iregge < 0).
class ReggePick {
  public:
    bool pick(Rng& rng, double ppbeam, double bslope, double tcutoff,
              double sminf, double rmtarget, double rmproduce,
              double rmrecoil, double& cosinethetacm);

  private:
    bool first_ = true;
    double dsigmadtnot_ = 1000.0, dsdtmax_ = 1000.0, dsdtmin_ = 0.0;
    double bslopelocal_ = 0.0;
};

// DOUBLEREGGEPICK: accept/reject the bottom-vertex t of the double-Regge
// mechanism against exp(2 b ln(s/s0) t).
bool doubleReggePick(Rng& rng, double ppbeam, double rmtarget, double tt,
                     double bslope2);

// SCULPTMASS: accept invariant mass rm23 against a shape: ipick=1
// exponential in (m - m_min) with slope bslope2; ipick=2 S-wave scattering
// length (bslope2) + effective range (bslope3) form.
bool sculptMass(Rng& rng, double rm23, double rmassmin, double bslope2,
                double bslope3, int ipick);

// PHASESPACEWEIGHT: keep events proportionally to the 2-body decay momentum
// (adaptive maximum, as in the Fortran).
class PhaseSpaceWeight {
  public:
    bool keep(Rng& rng, double m0, double m1, double m2, double pcmin);

  private:
    double psratiomax_ = 0.0;
};

// PHASESPACEWEIGHT3: sequential-2-body model of a 3-body final state.
class PhaseSpaceWeight3 {
  public:
    bool keep(Rng& rng, double m0, double m1, double m2, double pcminter,
              double pcmvaluein, double wvalue);

  private:
    bool first_ = true;
    double threeweightmax_ = 0.5, threeweightmin_ = 0.0;
};

// IHAVESPINHALF list from the Fortran.
bool haveSpinHalf(int lundId);

// POLARIZE: toy-model polarization axes for strongly produced particles.
void polarize(const FourVec& pv, const FourVec& p1, const FourVec& p2,
              int lundid1, int lundid2, Vec3& pol1, Vec3& pol2, double cthcm);

// WEAKDECAY: hyperon-style weak decay with asymmetry alpha = 0.642
// correlated with the parent polarization polv.
bool weakDecay(Rng& rng, const FourVec& p0, double m1, double m2,
               const Vec3& polv, FourVec& p1, FourVec& p2, Vec3& pol1,
               Vec3& pol2, double& cthcm, int lund);

// PRODUCT_STUFFIT: two-body decay of p0 with a pre-selected cos(theta_cm)
// (relative to the beam direction) and a random azimuth; the SECOND product
// heads along the chosen direction (recoil listed first in .def files).
bool productStuffit(Rng& rng, const FourVec& p0, double m1, double m2,
                    FourVec& p1, FourVec& p2, double& cthcm);

// PRODUCT_ROTATE: orient a CM-frame 3-body system so the composite p23
// points at costhetacm (random azimuth).  Vectors are modified in place;
// returns false when the final check fails (rare numerical issue).
bool productRotate(Rng& rng, const FourVec& p0, FourVec& p1, FourVec& p2,
                   FourVec& p3, FourVec& p23, double costhetacm);

// BOOST_TO_LAB for the product_rotate system.
void boostToLab(const FourVec& p0, FourVec& p1, FourVec& p2, FourVec& p3,
                FourVec& p23);

// GETT: Mandelstam t between pa and pb, plus tmin (parallel momenta) and
// tprime = t - tmin.
void gett(const FourVec& pa, const FourVec& pb, double& tt, double& ttmin,
          double& ttprime);

// INVARMASS: pt = p1 + p2 with the invariant mass in pt.m; returns the mass.
double invarMass(const FourVec& p1, const FourVec& p2, FourVec& pt);

// GSMEAR: crude GlueX instrumental smearing.  Photons (lund==1) get
// calorimeter energy smearing; everything else gets per-component momentum
// smearing with dp/p = dpoverp (+ small momentum-dependent term).
FourVec gsmear(Rng& rng, const FourVec& ppure, int lund, double dpoverp);

}  // namespace mcgen
