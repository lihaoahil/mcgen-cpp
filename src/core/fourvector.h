#pragma once
// Four-vector and kinematics utilities ported from MC_GEN's kinpack2.F and
// carrot.F (R. A. Schumacher, CMU).
//
// Convention preserved from the Fortran: a "four-vector" carries five
// components (px, py, pz, mass, energy).  The Fortran used single-precision
// REAL; this port uses double throughout (intentional deviation, documented
// in the project README).

#include <array>
#include <cmath>

#include "core/rng.h"

namespace mcgen {

struct FourVec {
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    double m  = 0.0;  // invariant mass
    double e  = 0.0;  // total energy

    double p2() const { return px * px + py * py + pz * pz; }          // PSQ
    double p() const { return std::sqrt(p2()); }                       // ABSP
    double etot() const { return std::sqrt(p2() + m * m); }            // ETOT
    void   updateEnergy() { e = etot(); }
};

using Vec3 = std::array<double, 3>;

inline double dot3(const FourVec& a, const FourVec& b) {              // DOT
    return a.px * b.px + a.py * b.py + a.pz * b.pz;
}

// EPCM: total momentum and invariant mass (stored in .m) of a two-particle
// system; .e holds the summed lab energy.
FourVec epcm(const FourVec& p1, const FourVec& p2);

// KIN12: decay of particle p0 into two bodies with masses m3, m4, with the
// first daughter emitted along direction dir3 (need not be normalized).
// Chooses the higher-momentum solution when two exist.  Returns false when
// the decay is kinematically forbidden (Fortran IERR=1).
bool kin12(const FourVec& p0, double m3, double m4, const Vec3& dir3,
           FourVec& p3, FourVec& p4);

// KIN22: 2-body reaction p1 + p2 -> p3 + p4 with p3 along dir3.
bool kin22(const FourVec& p1, const FourVec& p2, double m3, double m4,
           const Vec3& dir3, FourVec& p3, FourVec& p4);

// BOOST: boost p1 into the frame moving with velocity beta (units of c).
FourVec boost(const FourVec& p1, const Vec3& beta);

// ROTATE / UNROT: rotate a vector into / out of the frame in which dir3
// points along +z.
FourVec rotate(const FourVec& p1, const Vec3& dir3);
FourVec unrot(const FourVec& p1, const Vec3& dir3);

// CARROT: rotate the 3-vector part of vec by rotang radians about axis rot
// (clockwise looking along the axis).  No-op when the axis is null.
void carrot(FourVec& vec, const FourVec& rot, double rotang);

// CROSSPROD: c = a x b; returns the angle between a and b in radians.
double crossprod(const FourVec& a, const FourVec& b, FourVec& c);

// JACOB: boost along +z with speed bet; returns boosted vector, lab angle
// theta (degrees) and the Jacobian d(cos theta_cm)/d(cos theta_lab).
FourVec jacob(const FourVec& p1, double bet, double& theta, double& ajac);

// TRAVEL: straight-line propagation rout = rin + d * dir.
Vec3 travel(const Vec3& rin, const Vec3& dir, double d);

// KIN1R: isotropic two-body decay of p0 into masses m1, m2, boosted to the
// lab.  Outputs the daughters, the decay cosine w.r.t. the p0 direction
// (cthcm) and the rest-frame daughter momentum (pcm).
bool kin1r(Rng& rng, const FourVec& p0, double m1, double m2, FourVec& p1,
           FourVec& p2, double& cthcm, double& pcm);

// KIN2R: S-wave 2-body reaction p1 + p2 -> p3 + p4 at a random angle.
bool kin2r(Rng& rng, const FourVec& p1, const FourVec& p2, double m3,
           double m4, FourVec& p3, FourVec& p4, double& cthp3);

// KIN1RM: decay two parents (p0 -> p1+p2, p3 -> p4+p5) with the identical
// randomly chosen rest-frame direction, each boosted to the lab.
bool kin1rm(Rng& rng, const FourVec& p0, const FourVec& p3, double m1,
            double m2, double m4, double m5, FourVec& p1, FourVec& p2,
            FourVec& p4, FourVec& p5, double& cthcm);

}  // namespace mcgen
