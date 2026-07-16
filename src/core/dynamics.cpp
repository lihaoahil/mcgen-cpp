#include "core/dynamics.h"

#include <cmath>
#include <cstdio>

namespace mcgen {

namespace {
constexpr double kTwoPi = 2.0 * M_PI;

// Two-body CM energy of particle a in a system of mass w against b.
double cmEnergy(double w, double ma2minusmb2) {
    return (w * w - ma2minusmb2) / (2.0 * w);
}
}  // namespace

bool ReggePick::pick(Rng& rng, double ppbeam, double bslope, double tcutoff,
                     double sminf, double rmtarget, double rmproduce,
                     double rmrecoil, double& cosinethetacm) {
    const double snot = 1.0;  // no reliable estimate; default (see Fortran)
    if (first_ || bslope != bslopelocal_) {
        dsigmadtnot_ = 1000.0;
        dsdtmax_ = dsigmadtnot_;
        dsdtmin_ = dsigmadtnot_ * sminf;
        if (first_) {
            first_ = false;
            std::printf(" Regge exchange parameters:\n"
                        " t-slope (b)           = %8.3f GeV^-2\n"
                        " t-slope cutoff        = %8.3f GeV^2\n"
                        " Scale factor s_not    = %8.3f GeV^2\n",
                        bslope, tcutoff, snot);
        }
        bslopelocal_ = bslope;
    }

    const double rmg = 0.0;  // photon mass
    const double ww = std::sqrt(2.0 * ppbeam * rmtarget + rmtarget * rmtarget);
    const double egcm = cmEnergy(ww, rmtarget * rmtarget - rmg * rmg);
    const double eproducecm =
        cmEnergy(ww, rmrecoil * rmrecoil - rmproduce * rmproduce);
    const double pgcm = egcm;
    double pproducecm = eproducecm * eproducecm - rmproduce * rmproduce;
    if (pproducecm < 0.0 && pproducecm > -2.0e-5)
        pproducecm = 1.0e-6;  // round-off fix: place right at threshold
    if (pproducecm >= 1.0e-8) {
        pproducecm = std::sqrt(pproducecm);
    } else {
        std::printf(" CM momentum incorrect: %g %g %g %g\n", ppbeam, ww,
                    rmtarget, pproducecm);
        return false;
    }
    const double m2sum = rmproduce * rmproduce + rmg * rmg;
    const double tmin = m2sum - 2.0 * (egcm * eproducecm - pgcm * pproducecm);
    const double tmax = m2sum - 2.0 * (egcm * eproducecm + pgcm * pproducecm);
    const double ss = ww * ww;

    double tt = 0.0;
    int iregge = 0;
    for (;;) {
        ++iregge;
        const double tprime = rng.uniform() * (tmax - tmin);
        tt = tmin + tprime;
        double dsdt;
        if (std::abs(tprime) > std::abs(tcutoff)) {
            dsdt = dsigmadtnot_ *
                   std::exp(2.0 * bslope * std::log(ss / snot) *
                            (tprime + tcutoff));
        } else {
            // phenomenological linear rise at low t'
            dsdt = (tcutoff != 0.0)
                       ? dsigmadtnot_ * std::abs(tprime / tcutoff)
                       : dsigmadtnot_;
        }
        dsdt = std::max(dsdt, dsdtmin_);
        if (dsdt > dsdtmax_) {
            std::printf(" Raising dsdtmax: %g %g\n", dsdt, dsdtmax_);
            dsdtmax_ = dsdt;
        }
        if (rng.uniform() * dsdtmax_ <= dsdt) break;
        if (iregge >= 2000 && iregge % 2000 == 0)
            std::printf(" Inefficient t-channel generation: %d %g %g\n",
                        iregge, tmin, tmax);
    }
    cosinethetacm =
        (tt - m2sum + 2.0 * (egcm * eproducecm)) / (2.0 * pgcm * pproducecm);
    if (std::isnan(cosinethetacm)) {
        std::printf(" REGGEPICK produced NaN cosine theta\n");
        return false;
    }
    if (cosinethetacm > 1.0) cosinethetacm = 1.0;    // yes, this can happen
    if (cosinethetacm < -1.0) cosinethetacm = -1.0;
    return true;
}

bool doubleReggePick(Rng& rng, double ppbeam, double rmtarget, double tt,
                     double bslope2) {
    const double scale = 1000.0;
    const double ss = 2.0 * ppbeam * rmtarget + rmtarget * rmtarget;
    const double snot = 1.0;
    const double probdensity =
        scale * std::exp(2.0 * bslope2 * std::log(ss / snot) * tt);
    return scale * rng.uniform() < probdensity;
}

bool sculptMass(Rng& rng, double rm23, double rmassmin, double bslope2,
                double bslope3, int ipick) {
    const double hbarc = 197.0;  // MeV fm
    const double xx = rm23 - rmassmin;
    double yy, ymax;
    if (ipick == 1) {  // exponential shape
        yy = std::exp(-xx / bslope2);
        ymax = 1.0;
    } else if (ipick == 2) {  // scattering length + effective range
        if (bslope2 == 0.0) return false;
        ymax = bslope2 * bslope2;
        const double qq = std::sqrt((0.5 * rm23) * (0.5 * rm23) -
                                    (0.5 * rmassmin) * (0.5 * rmassmin));
        const double wk = qq / hbarc * 1000.0;  // wavenumber in 1/fm
        const double d = 1.0 - 0.5 * bslope2 * bslope3 * wk * wk;
        yy = ymax / (bslope2 * wk * wk + d * d);
    } else {
        std::printf(" SCULPTMASS: no such option %d\n", ipick);
        std::exit(1);
    }
    return rng.uniform() * ymax < yy;
}

bool PhaseSpaceWeight::keep(Rng& rng, double m0, double m1, double m2,
                            double pcmin) {
    const double a = m0 * m0 - (m1 * m1 + m2 * m2);
    double pcmvalue = a * a - (2.0 * m1 * m2) * (2.0 * m1 * m2);
    pcmvalue = std::sqrt(pcmvalue) / (2.0 * m0);
    const double psratio = pcmvalue / pcmin;
    if (psratio > psratiomax_) psratiomax_ = psratio;
    return rng.uniform() * psratiomax_ <= psratio;
}

bool PhaseSpaceWeight3::keep(Rng& rng, double m0, double m1, double m2,
                             double pcminter, double pcmvaluein,
                             double wvalue) {
    if (first_) {
        first_ = false;
        std::printf(" First call to PHASESPACEWEIGHT3\n");
    }
    const double a = m0 * m0 - (m1 * m1 + m2 * m2);
    double pcmvalue = a * a - (2.0 * m1 * m2) * (2.0 * m1 * m2);
    pcmvalue = std::sqrt(pcmvalue) / (2.0 * m0);
    double threeweight =
        (pcminter * pcmvalue) / (pcmvaluein * wvalue * wvalue) * 10.0;
    if (threeweight > threeweightmax_) {
        threeweightmax_ = threeweight;
        std::printf(" PS3> New maximum phase space %g\n", threeweight);
    }
    if (threeweight < threeweightmin_) {
        threeweightmin_ = threeweight;
        std::printf(" PS3> New minimum phase space %g\n", threeweight);
    }
    return rng.uniform() * threeweightmax_ <= threeweight;
}

bool haveSpinHalf(int lundId) {
    static constexpr int kList[] = {999, 57, 41, 42, 43, 44, 45,
                                    46,  47, 49, 94, 7,  2,  9};
    const int it = std::abs(lundId);
    for (int v : kList)
        if (it == v) return true;
    return false;
}

void polarize(const FourVec& pv, const FourVec& p1, const FourVec& p2,
              int lundid1, int lundid2, Vec3& pol1, Vec3& pol2, double cthcm) {
    FourVec khat;
    int iflag;
    if (haveSpinHalf(lundid1) && haveSpinHalf(lundid2)) {
        crossprod(p1, p2, khat);
        iflag = 2;  // polarize both along their mutual normal
    } else if (haveSpinHalf(lundid1)) {
        crossprod(pv, p1, khat);
        iflag = 11;
    } else if (haveSpinHalf(lundid2)) {
        crossprod(pv, p2, khat);
        iflag = 12;
    } else {
        pol1 = {0, 0, 0};
        pol2 = {0, 0, 0};
        return;
    }
    double rkhatnorm = khat.p();
    if (rkhatnorm == 0.0) rkhatnorm = 1.e36;
    if (cthcm > 1.0) {
        std::printf(" Polarize angle bad %g\n", cthcm);
        cthcm = 1.0;
    }
    // Toy model mimicking kaon photoproduction.
    const double fact = std::sqrt(1.0 - cthcm * cthcm) * (0.3 + 1.0 * cthcm);
    if (iflag == 2) {
        const double f = 1.0 / rkhatnorm;
        pol1 = {khat.px * f, khat.py * f, khat.pz * f};
        pol2 = {-pol1[0], -pol1[1], -pol1[2]};
    } else if (iflag == 11) {
        const double f = fact / rkhatnorm;
        pol1 = {khat.px * f, khat.py * f, khat.pz * f};
        pol2 = {0, 0, 0};
    } else {
        const double f = fact / rkhatnorm;
        pol2 = {khat.px * f, khat.py * f, khat.pz * f};
        pol1 = {0, 0, 0};
    }
}

bool weakDecay(Rng& rng, const FourVec& p0, double m1, double m2,
               const Vec3& polv, FourVec& p1, FourVec& p2, Vec3& pol1,
               Vec3& pol2, double& cthcm, int lund) {
    const double alphalam = 0.642;
    double alpha;
    if (lund == 57 || lund == -57) {
        // anti-Lambda: opposite asymmetry, but POLV is already flipped
        alpha = alphalam;
    } else {
        std::printf(" Decay asymmetry is not specified %d\n", lund);
        alpha = 1.0;
    }

    double polnorm = std::sqrt(polv[0] * polv[0] + polv[1] * polv[1] +
                               polv[2] * polv[2]);
    if (polnorm > 1.0) polnorm = 1.0;
    const double a = alpha * polnorm;

    const double y = rng.uniform();
    // Linearly weighted cosine of the polar angle w.r.t. polv.
    double x = (a != 0.0)
                   ? (std::sqrt((a - 1.0) * (a - 1.0) + 4.0 * a * y) - 1.0) / a
                   : 2.0 * y - 1.0;
    if (x > 1.0) x = 1.0;
    if (x < -1.0) x = -1.0;
    const double theta = std::acos(x);
    const double phik = kTwoPi * rng.uniform();

    // Decay direction: start along POLV, rotate by theta about (polv x p0),
    // then by phik about polv.
    FourVec dtmp;
    dtmp.px = polv[0];
    dtmp.py = polv[1];
    dtmp.pz = polv[2];
    FourVec polv4 = dtmp;
    FourVec rot;
    crossprod(polv4, p0, rot);
    carrot(dtmp, rot, theta);
    carrot(dtmp, polv4, phik);

    cthcm = dot3(dtmp, p0) / (dtmp.p() * p0.p());

    const double e0 = p0.e;
    const Vec3 beta{-p0.px / e0, -p0.py / e0, -p0.pz / e0};
    FourVec p0rest;
    p0rest.m = p0.m;
    p0rest.updateEnergy();
    FourVec p1tmp, p2tmp;
    // Particle 2 heads along DTMP (Fortran passes p1tmp first but the decay
    // direction argument applies to the first product = p1tmp there; the
    // call is kin12(p0tmp, p1tmp, p2tmp, dtmp) so p1 goes along dtmp).
    if (!kin12(p0rest, m1, m2, {dtmp.px, dtmp.py, dtmp.pz}, p1tmp, p2tmp))
        return false;
    p1 = boost(p1tmp, beta);
    p2 = boost(p2tmp, beta);
    pol1 = {0, 0, 0};
    pol2 = {0, 0, 0};
    return true;
}

bool productStuffit(Rng& rng, const FourVec& p0, double m1, double m2,
                    FourVec& p1, FourVec& p2, double& cthcm) {
    const double theta = std::acos(cthcm);
    const double phik = kTwoPi * rng.uniform();

    FourVec ex;
    ex.px = 1.0;
    FourVec dtmp = p0;
    FourVec rot;
    crossprod(p0, ex, rot);
    carrot(dtmp, rot, theta);   // polar rotation
    carrot(dtmp, p0, phik);     // azimuthal rotation

    cthcm = dot3(dtmp, p0) / (dtmp.p() * p0.p());
    if (std::abs(cthcm - std::cos(theta)) > 0.001)
        std::printf(" Deep trouble: %g %g\n", cthcm, std::cos(theta));

    const Vec3 beta{-p0.px / p0.e, -p0.py / p0.e, -p0.pz / p0.e};
    FourVec p0rest;
    p0rest.m = p0.m;
    p0rest.updateEnergy();
    FourVec p1tmp, p2tmp;
    // iproduct=2: the SECOND product heads along dtmp (recoil listed first).
    if (!kin12(p0rest, m2, m1, {dtmp.px, dtmp.py, dtmp.pz}, p2tmp, p1tmp))
        return false;
    p1 = boost(p1tmp, beta);
    p2 = boost(p2tmp, beta);
    return true;
}

bool productRotate(Rng& rng, const FourVec& p0, FourVec& p1, FourVec& p2,
                   FourVec& p3, FourVec& p23, double costhetacm) {
    const double theta = std::acos(costhetacm);
    const double phik = kTwoPi * rng.uniform();

    FourVec ex, ze;
    ex.px = 1.0;
    ze.pz = 1.0;

    // Rotate the whole CM system so p23 points along the beam axis.
    FourVec rot;
    const double angle = crossprod(p23, ze, rot);
    carrot(p23, rot, angle);
    carrot(p1, rot, angle);
    carrot(p2, rot, angle);
    carrot(p3, rot, angle);

    // Tilt to the desired polar angle (about X), then scramble the azimuth
    // (about the beam = p0 direction).
    carrot(p23, ex, theta);
    carrot(p1, ex, theta);
    carrot(p2, ex, theta);
    carrot(p3, ex, theta);
    carrot(p23, p0, phik);
    carrot(p1, p0, phik);
    carrot(p2, p0, phik);
    carrot(p3, p0, phik);

    const double cthcm = dot3(p23, ze) / (p23.p() * 1.0);
    return std::abs(cthcm - std::cos(theta)) <= 0.001;
}

void boostToLab(const FourVec& p0, FourVec& p1, FourVec& p2, FourVec& p3,
                FourVec& p23) {
    const Vec3 beta{-p0.px / p0.e, -p0.py / p0.e, -p0.pz / p0.e};
    p1 = boost(p1, beta);
    p2 = boost(p2, beta);
    p3 = boost(p3, beta);
    p23 = boost(p23, beta);
}

void gett(const FourVec& pa, const FourVec& pb, double& tt, double& ttmin,
          double& ttprime) {
    const double dpx = pb.px - pa.px, dpy = pb.py - pa.py, dpz = pb.pz - pa.pz;
    tt = (pb.e - pa.e) * (pb.e - pa.e) - (dpx * dpx + dpy * dpy + dpz * dpz);

    const double pmagp = pb.p();
    FourVec p0;
    p0.pz = (pa.pz > 0.0) ? pmagp : -pmagp;
    p0.m = pb.m;
    p0.e = pb.e;
    const double d0z = p0.pz - pa.pz;
    ttmin = (p0.e - pa.e) * (p0.e - pa.e) -
            (pa.px * pa.px + pa.py * pa.py + d0z * d0z);
    ttprime = tt - ttmin;
}

double invarMass(const FourVec& p1, const FourVec& p2, FourVec& pt) {
    pt.px = p1.px + p2.px;
    pt.py = p1.py + p2.py;
    pt.pz = p1.pz + p2.pz;
    pt.e = p1.e + p2.e;
    const double m2 = pt.e * pt.e - pt.p2();
    if (m2 <= 0.0) {
        std::printf(" INVARMASS Ack!\n");
        std::exit(1);
    }
    pt.m = std::sqrt(m2);
    return pt.m;
}

FourVec gsmear(Rng& rng, const FourVec& ppure, int lund, double dpoverp) {
    FourVec pspect = ppure;
    const double pmom = ppure.p();
    const double evalue = std::sqrt(pmom * pmom + ppure.m * ppure.m);
    if (lund != 1) {  // "charged" particle: momentum + angle smearing
        if (dpoverp != 0.0) {
            auto smear1 = [&](double comp) {
                const double dpop =
                    std::sqrt(dpoverp * dpoverp + (0.005 * comp) * (0.005 * comp));
                return comp * (1.0 + dpop * rng.normal());
            };
            pspect.px = smear1(ppure.px);
            pspect.py = smear1(ppure.py);
            pspect.pz = smear1(ppure.pz);
        }
    } else {  // photon: energy smearing only
        double theta_polar = 1000.0;
        if (pmom != 0.0) {
            double arg = ppure.pz / pmom;
            arg = std::clamp(arg, -1.0, 1.0);
            theta_polar = std::acos(arg) * 180.0 / M_PI;
        }
        const double deoe = (theta_polar < 11.0)
                                ? 0.05 / std::sqrt(evalue)    // FCAL
                                : 0.135 / std::sqrt(evalue);  // BCAL
        const double ratio = 1.0 + deoe * rng.normal();
        pspect.px = ppure.px * ratio;
        pspect.py = ppure.py * ratio;
        pspect.pz = ppure.pz * ratio;
    }
    pspect.m = ppure.m;
    pspect.updateEnergy();
    return pspect;
}

}  // namespace mcgen
