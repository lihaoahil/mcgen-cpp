#include "core/fourvector.h"

namespace mcgen {

FourVec epcm(const FourVec& p1, const FourVec& p2) {
    FourVec pcm;
    pcm.px = p1.px + p2.px;
    pcm.py = p1.py + p2.py;
    pcm.pz = p1.pz + p2.pz;
    const double e = p1.etot() + p2.etot();
    pcm.m = std::sqrt(e * e - pcm.p2());  // invariant mass W
    pcm.e = e;
    return pcm;
}

bool kin12(const FourVec& p0, double m3, double m4, const Vec3& dir3,
           FourVec& p3, FourVec& p4) {
    const double ap0s = p0.p2();
    const double adir = std::sqrt(dir3[0] * dir3[0] + dir3[1] * dir3[1] +
                                  dir3[2] * dir3[2]);
    double calph = 1.0;
    if (ap0s != 0.0) {
        calph = (p0.px * dir3[0] + p0.py * dir3[1] + p0.pz * dir3[2]) /
                std::sqrt(ap0s) / adir;
    }
    const double e0s = p0.m * p0.m + ap0s;
    const double d = e0s + m3 * m3 - m4 * m4 - ap0s;
    const double a = 4.0 * (ap0s * calph * calph - e0s);
    const double b = 4.0 * d * std::sqrt(ap0s) * calph;
    const double c = d * d - 4.0 * e0s * m3 * m3;
    const double disc = b * b - 4.0 * a * c;

    p3.m = m3;
    p4.m = m4;
    if (disc < 0.0) return false;

    // Higher-momentum solution, as in the Fortran.
    const double ap3 = (-b - std::sqrt(disc)) / (2.0 * a);
    p3.px = ap3 * dir3[0] / adir;
    p3.py = ap3 * dir3[1] / adir;
    p3.pz = ap3 * dir3[2] / adir;
    p4.px = p0.px - p3.px;
    p4.py = p0.py - p3.py;
    p4.pz = p0.pz - p3.pz;
    p3.updateEnergy();
    p4.updateEnergy();
    return true;
}

bool kin22(const FourVec& p1, const FourVec& p2, double m3, double m4,
           const Vec3& dir3, FourVec& p3, FourVec& p4) {
    return kin12(epcm(p1, p2), m3, m4, dir3, p3, p4);
}

namespace {
// Shared by ROTATE and UNROT: the two rotation angles that carry dir3 to +z.
struct RotAngles {
    double sina1, cosa1, sina2, cosa2;
};

RotAngles rotAngles(const Vec3& dir3) {
    RotAngles a;
    double r = std::sqrt(dir3[0] * dir3[0] + dir3[1] * dir3[1]);
    if (r == 0.0) {
        a.sina1 = 0.0;
        a.cosa1 = 1.0;
    } else {
        a.sina1 = -dir3[0] / r;
        a.cosa1 = dir3[1] / r;
    }
    const double ytmp = -a.sina1 * dir3[0] + a.cosa1 * dir3[1];
    r = std::sqrt(ytmp * ytmp + dir3[2] * dir3[2]);
    if (r == 0.0) {
        a.sina2 = 0.0;
        a.cosa2 = 1.0;
    } else {
        a.sina2 = -ytmp / r;
        a.cosa2 = dir3[2] / r;
    }
    return a;
}
}  // namespace

FourVec rotate(const FourVec& p1, const Vec3& dir3) {
    const RotAngles a = rotAngles(dir3);
    FourVec p2 = p1;
    p2.px = a.cosa1 * p1.px + a.sina1 * p1.py;
    const double ytmp = -a.sina1 * p1.px + a.cosa1 * p1.py;
    p2.py = a.cosa2 * ytmp + a.sina2 * p1.pz;
    p2.pz = -a.sina2 * ytmp + a.cosa2 * p1.pz;
    return p2;
}

FourVec unrot(const FourVec& p1, const Vec3& dir3) {
    const RotAngles a = rotAngles(dir3);
    FourVec p2 = p1;
    p2.pz = a.sina2 * p1.py + a.cosa2 * p1.pz;
    const double ytmp = a.cosa2 * p1.py - a.sina2 * p1.pz;
    p2.py = a.sina1 * p1.px + a.cosa1 * ytmp;
    p2.px = a.cosa1 * p1.px - a.sina1 * ytmp;
    return p2;
}

FourVec boost(const FourVec& p1, const Vec3& beta) {
    const double bet = std::sqrt(beta[0] * beta[0] + beta[1] * beta[1] +
                                 beta[2] * beta[2]);
    const double gam = 1.0 / std::sqrt(1.0 - bet * bet);
    const double e1 = p1.e;

    FourVec ptmp = rotate(p1, beta);
    ptmp.pz = gam * (ptmp.pz - bet * e1);
    FourVec p2 = unrot(ptmp, beta);
    p2.m = p1.m;
    p2.updateEnergy();
    return p2;
}

void carrot(FourVec& vec, const FourVec& rot, double rotang) {
    const double rr = rot.p();
    if (rr == 0.0) return;
    const Vec3 rhat{rot.px / rr, rot.py / rr, rot.pz / rr};

    const double vecdot = vec.px * rhat[0] + vec.py * rhat[1] + vec.pz * rhat[2];
    const Vec3 vecpara{vecdot * rhat[0], vecdot * rhat[1], vecdot * rhat[2]};
    FourVec vecperp;
    vecperp.px = vec.px - vecpara[0];
    vecperp.py = vec.py - vecpara[1];
    vecperp.pz = vec.pz - vecpara[2];

    FourVec rhat4;
    rhat4.px = rhat[0];
    rhat4.py = rhat[1];
    rhat4.pz = rhat[2];
    FourVec vcross;
    crossprod(vecperp, rhat4, vcross);

    const double sinrot = std::sin(rotang);
    const double cosrot = std::cos(rotang);
    vec.px = vecperp.px * cosrot - vcross.px * sinrot + vecpara[0];
    vec.py = vecperp.py * cosrot - vcross.py * sinrot + vecpara[1];
    vec.pz = vecperp.pz * cosrot - vcross.pz * sinrot + vecpara[2];
}

double crossprod(const FourVec& a, const FourVec& b, FourVec& c) {
    c.px = b.pz * a.py - b.py * a.pz;
    c.py = -(b.pz * a.px - b.px * a.pz);
    c.pz = b.py * a.px - b.px * a.py;
    c.m = 0.0;
    c.e = 0.0;
    const double aa = a.p2();
    const double bb = b.p2();
    const double dd = dot3(a, b);
    return std::acos(dd / std::sqrt(aa * bb));
}

FourVec jacob(const FourVec& p1, double bet, double& theta, double& ajac) {
    const double gam = 1.0 / std::sqrt(1.0 - bet * bet);
    const double e1 = p1.etot();
    FourVec p2 = p1;
    p2.pz = gam * p1.pz + bet * gam * e1;
    p2.updateEnergy();
    const double p2abs = p2.p();
    const double bet2 = p2abs / p2.etot();
    const double cos2 = p2.pz / p2abs;
    theta = 180.0 / M_PI * std::acos(cos2);
    if (bet == 0.0) {
        ajac = 1.0;
    } else {
        const double gamsq = 1.0 / (1.0 - bet * bet);
        ajac = (1.0 - bet / bet2 * cos2) /
               (gamsq * std::pow((bet / bet2 - cos2) * (bet / bet2 - cos2) +
                                     (1.0 - cos2 * cos2) / gamsq,
                                 1.5));
    }
    return p2;
}

Vec3 travel(const Vec3& rin, const Vec3& dir, double d) {
    return {rin[0] + d * dir[0], rin[1] + d * dir[1], rin[2] + d * dir[2]};
}

namespace {
Vec3 randomDirection(Rng& rng, double& cthcm) {
    cthcm = 1.0 - 2.0 * rng.uniform();
    const double phi = 2.0 * M_PI * rng.uniform();
    const double tmp = 1.0 - cthcm * cthcm;
    const double sthcm = tmp > 0.0 ? std::sqrt(tmp) : 0.0;
    return {sthcm * std::cos(phi), sthcm * std::sin(phi), cthcm};
}
}  // namespace

bool kin1r(Rng& rng, const FourVec& p0, double m1, double m2, FourVec& p1,
           FourVec& p2, double& cthcm, double& pcm) {
    const Vec3 dtmp = randomDirection(rng, cthcm);
    const double e0 = p0.etot();
    const Vec3 beta{-p0.px / e0, -p0.py / e0, -p0.pz / e0};

    FourVec p0rest;
    p0rest.m = p0.m;
    p0rest.updateEnergy();
    FourVec p1tmp, p2tmp;
    if (!kin12(p0rest, m1, m2, dtmp, p1tmp, p2tmp)) return false;
    p1 = boost(p1tmp, beta);
    p2 = boost(p2tmp, beta);
    pcm = p1tmp.p();
    return true;
}

bool kin2r(Rng& rng, const FourVec& p1, const FourVec& p2, double m3,
           double m4, FourVec& p3, FourVec& p4, double& cthp3) {
    double pcm = 0.0;
    return kin1r(rng, epcm(p1, p2), m3, m4, p3, p4, cthp3, pcm);
}

bool kin1rm(Rng& rng, const FourVec& p0, const FourVec& p3, double m1,
            double m2, double m4, double m5, FourVec& p1, FourVec& p2,
            FourVec& p4, FourVec& p5, double& cthcm) {
    const Vec3 dtmp = randomDirection(rng, cthcm);
    const double e0 = p0.etot();
    const double e3 = p3.etot();
    const Vec3 beta0{-p0.px / e0, -p0.py / e0, -p0.pz / e0};
    const Vec3 beta3{-p3.px / e3, -p3.py / e3, -p3.pz / e3};

    FourVec p0rest, p3rest;
    p0rest.m = p0.m;
    p0rest.updateEnergy();
    p3rest.m = p3.m;
    p3rest.updateEnergy();

    FourVec p1tmp, p2tmp, p4tmp, p5tmp;
    if (!kin12(p0rest, m1, m2, dtmp, p1tmp, p2tmp)) return false;
    if (!kin12(p3rest, m4, m5, dtmp, p4tmp, p5tmp)) return false;
    p1 = boost(p1tmp, beta0);
    p2 = boost(p2tmp, beta0);
    p4 = boost(p4tmp, beta3);
    p5 = boost(p5tmp, beta3);
    return true;
}

}  // namespace mcgen
