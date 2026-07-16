#include "core/targetgeom.h"

#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace mcgen {

namespace {
// LINECIRCLE: intersection of the (x0,y0)->(x1,y1) line with a circle.
void linecircle(double x0, double y0, double x1, double y1, double radius,
                double& xx, double& yy) {
    const double theta = std::atan2(y1 - y0, x1 - x0);
    const double cost = std::cos(theta);
    const double sint = std::sin(theta);
    const double b = 2.0 * (y0 * sint + x0 * cost);
    const double c = x0 * x0 + y0 * y0 - radius * radius;
    const double rad = std::sqrt(b * b - 4.0 * c);
    const double tplus = 0.5 * (-b + rad);
    xx = x0 + tplus * cost;
    yy = y0 + tplus * sint;
}

// Classify a point: 0 well inside the inner radius, 1 in the wall shell,
// 2 outside.
int classify(const Vec3& r, double radius1, double radius2, double zmin,
             double zmax) {
    constexpr double eps = 0.001;
    const double rr = std::sqrt(r[0] * r[0] + r[1] * r[1]);
    if (rr < radius2 && r[2] >= zmin && r[2] <= zmax) {
        if (rr < radius1 - eps && r[2] >= zmin + eps && r[2] <= zmax - eps)
            return 0;
        return 1;
    }
    return 2;
}
}  // namespace

int checkout(const Vec3& r0, Vec3& r1, const Vec3& direc,
             const TargetGeom& g) {
    if (g.icode == 0) return 0;
    if (g.icode != 2)
        throw std::runtime_error(
            "checkout: only cylindrical target geometry is supported "
            "(as in the Fortran)");
    const double radius1 = g.size1 * 0.5;
    const double radius2 = g.size2 * 0.5;
    const double zmin = g.size3, zmax = g.size4;

    const int istart = classify(r0, radius1, radius2, zmin, zmax);
    const int iend = classify(r1, radius1, radius2, zmin, zmax);

    int iflag;
    double clipRadius;
    if (istart == 0 && iend == 0) return 0;
    if (istart == 0 && (iend == 1 || iend == 2)) {
        iflag = 1;
        clipRadius = radius1;
    } else if (istart == 1 && iend == 2) {
        iflag = 2;
        clipRadius = radius2;
    } else if (istart == 1 && iend == 1) {
        return 0;
    } else {
        return 3;  // outside or backward-traveling; ignore
    }

    double xx, yy;
    linecircle(r0[0], r0[1], r1[0], r1[1], clipRadius, xx, yy);
    double zz = r0[2];
    if (direc[0] != 0.0) zz = r0[2] + direc[2] / direc[0] * (xx - r0[0]);
    r1 = {xx, yy, zz};
    if (r1[2] > zmax) {
        const double ratio = (zmax - r0[2]) / (r1[2] - r0[2]);
        r1[2] = zmax;
        r1[0] = r0[0] + ratio * (r1[0] - r0[0]);
        r1[1] = r0[1] + ratio * (r1[1] - r0[1]);
    } else if (r1[2] < zmin) {
        const double ratio = (zmin - r0[2]) / (r1[2] - r0[2]);
        r1[2] = zmin;
        r1[0] = r0[0] + ratio * (r1[0] - r0[0]);
        r1[1] = r0[1] + ratio * (r1[1] - r0[1]);
    }
    return iflag;
}

Vec3 getVertex(Rng& rng, const TargetGeom& g) {
    if (g.jcode != 1) return {0.0, 0.0, 0.0};
    for (int ibad = 1; ibad <= 10; ++ibad) {
        const double x = rng.normal() * g.xsigma;
        const double y = rng.normal() * g.ysigma;
        const double z = g.size3 + (g.size4 - g.size3) * rng.uniform();
        Vec3 r0{x, y, z};
        Vec3 r1 = r0;
        Vec3 direc{-1.0, -1.0, -1.0};
        const double rmag = std::sqrt(x * x + y * y + z * z);
        if (rmag > 0.0) direc = {-x / rmag, -y / rmag, -z / rmag};
        if (checkout(r0, r1, direc, g) < 1) return r0;
        std::printf(" Vertex outside target %d\n", ibad);
    }
    throw std::runtime_error("getvertex: beam spot is really bad");
}

Vec3 propagate(Rng& rng, const FourVec& pv, const Vec3& r0, double ctau) {
    const double pvmag = pv.p();
    const double et = pv.etot();
    if (et <= 0.0) {
        std::printf(" Bogus particle in propagate\n");
        return r0;
    }
    Vec3 direc{0.0, 0.0, 0.0};
    if (pvmag > 0.0)
        direc = {pv.px / pvmag, pv.py / pvmag, pv.pz / pvmag};

    // 1-u keeps the argument in (0,1]; same distribution as the Fortran.
    const double tdecay = -ctau * std::log(1.0 - rng.uniform());
    const double beta = pvmag / et;
    const double gamma =
        (beta < 1.0) ? std::sqrt(1.0 / (1.0 - beta * beta)) : 100000.0;
    const double dstep = beta * gamma * tdecay;  // lab flight distance
    return travel(r0, direc, dstep);
}

}  // namespace mcgen
