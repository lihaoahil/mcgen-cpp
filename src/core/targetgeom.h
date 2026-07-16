#pragma once
// Primary-vertex selection and particle propagation, ported from
// getvertex (mc_gen.F) and propagate/checkout/linecircle (mc_propagate.F).

#include "core/fourvector.h"
#include "core/rng.h"

namespace mcgen {

struct TargetGeom {
    int icode = 0;      // {0,1,2} = {none, rectangular, cylindrical}
    double size1 = 0, size2 = 0, size3 = 0, size4 = 0;  // cm
    int jcode = 0;      // {0,1} = {point vertex, diffuse vertex}
    double xsigma = 0, ysigma = 0;  // cm
};

// GETVERTEX: pick the primary vertex; Gaussian in x/y, uniform in z within
// [size3, size4], rejected if outside the (cylindrical) target.
Vec3 getVertex(Rng& rng, const TargetGeom& g);

// CHECKOUT: clip the segment r0->r1 at the target boundary (cylindrical
// geometry only, as in the Fortran).  Returns the Fortran IFLAG:
// 0 inside, 1/2 crossing (r1 moved onto the boundary), 3 outside.
int checkout(const Vec3& r0, Vec3& r1, const Vec3& direc,
             const TargetGeom& g);

// PROPAGATE: straight-line flight over an exponentially distributed decay
// length (no energy loss, as in the current Fortran).  Returns the decay
// point r1.
Vec3 propagate(Rng& rng, const FourVec& pv, const Vec3& r0, double ctau);

}  // namespace mcgen
