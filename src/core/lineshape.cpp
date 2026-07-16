#include "core/lineshape.h"

#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace mcgen {

double bwran(Rng& rng) {
    // Inverse CDF of the Cauchy density truncated to [-30, 30].
    static const double lo = std::atan(-30.0);
    static const double hi = std::atan(30.0);
    return std::tan(lo + (hi - lo) * rng.uniform());
}

RelativisticBWSampler::RelativisticBWSampler(double rmnot, double gammanot,
                                             double rm1, double rm2,
                                             int lvalue)
    : rm1_(rm1), rm2_(rm2) {
    std::printf(" Initializing Relativistic Breit-Wigner %g %g\n", rm1, rm2);
    std::printf(" Orbital angular momentum is %d\n", lvalue);

    rmmin_ = rmnot - 4.0 * gammanot;
    const double rmmax = rmnot + 20.0 * gammanot;
    drm_ = (rmmax - rmmin_) / kNbin;

    auto rad = [&](double m) {
        const double a = m * m - (rm1 * rm1 + rm2 * rm2);
        return a * a - (2.0 * rm1 * rm2) * (2.0 * rm1 * rm2);
    };

    double pcmnot = 1.0;
    if (rad(rmnot) >= 0.0) {
        pcmnot = std::sqrt(rad(rmnot)) / (2.0 * rmnot);
    } else {
        std::printf(" Central mass value is below threshold %g %g %g\n"
                    " Falling back to a NON-relativistic Breit-Wigner\n",
                    rmnot, rm1, rm2);
        ignoreMasses_ = true;
    }

    cdf_.assign(kNbin + 1, 0.0);
    for (int i = 1; i <= kNbin; ++i) {
        const double rm = rmmin_ + (i - 1) * drm_;
        double bwvalue = 0.0;
        const double r = rad(rm);
        if (r >= 0.0 && !ignoreMasses_ && rm > 0.0) {
            // NB: the Fortran divides by 2*rmnot (not 2*rm) here; preserved.
            const double pcm = std::sqrt(r) / (2.0 * rmnot);
            const double gamma =
                gammanot * std::pow(pcm / pcmnot, 2.0 * lvalue + 1.0);
            const double d1 = rm * rm - rmnot * rmnot;
            bwvalue = 2.0 * rm * rmnot * gamma /
                      (d1 * d1 + (rmnot * gamma) * (rmnot * gamma));
        } else if (ignoreMasses_) {
            // Plain Lorentzian when the pole sits below threshold.
            const double d1 = rm * rm - rmnot * rmnot;
            bwvalue = 2.0 * rm * rmnot * gammanot /
                      (d1 * d1 + (rmnot * gammanot) * (rmnot * gammanot));
        }
        cdf_[static_cast<size_t>(i)] = cdf_[static_cast<size_t>(i) - 1] + bwvalue;
    }
    const double fact = cdf_[kNbin];
    if (fact <= 0.0)
        throw std::runtime_error(
            "RelativisticBWSampler: lineshape integrates to zero");
    for (auto& v : cdf_) v /= fact;
}

double RelativisticBWSampler::sample(Rng& rng) const {
    for (int attempt = 0; attempt < 20; ++attempt) {
        const double test = rng.uniform();
        double bwmass = 0.0;
        bool found = false;
        for (int i = 1; i <= kNbin; ++i) {
            if (cdf_[static_cast<size_t>(i)] >= test) {
                bwmass = rmmin_ + (i - 1) * drm_;
                // Blur within one bin to smooth binning artifacts.
                bwmass += drm_ * (rng.uniform() - 0.5);
                found = true;
                break;
            }
        }
        if (!found)
            throw std::runtime_error("RelativisticBWSampler: CDF lookup failed");
        if (bwmass >= rm1_ + rm2_ || ignoreMasses_) return bwmass;
    }
    throw std::runtime_error(
        "RelativisticBWSampler: below-threshold mismatch after 20 tries");
}

RelativisticLineshape::Result RelativisticLineshape::sample(
    Rng& rng, double rmnot, double gammanot, double rm1, double rm2,
    int lvalue, double rm3, double ww) {
    if (first_) {
        first_ = false;
        std::printf(" First call to Relativistic Breit-Wigner with phase "
                    "space: m0=%g Gamma=%g m1=%g m2=%g m3=%g W=%g L=%d\n",
                    rmnot, gammanot, rm1, rm2, rm3, ww, lvalue);
    }

    auto twoBodyP = [](double w, double ma, double mb) {
        const double a = w * w - (ma * ma + mb * mb);
        const double rad = a * a - (2.0 * ma * mb) * (2.0 * ma * mb);
        return rad >= 0.0 ? std::sqrt(rad) / (2.0 * w) : -1.0;
    };

    const double rmmin = std::max(rmnot - 4.0 * gammanot, rm1 + rm2);
    const double rmmax = std::min(ww - rm3, rmnot + 20.0 * gammanot);
    const double drm = rmmax - rmmin;
    if (drm < 0.0) return {rmnot, 1};

    const double pcmnot = twoBodyP(ww, rmmin, rm3);
    if (pcmnot < 0.0) return {rmnot, 4};
    const double qcmnot = twoBodyP(rmnot, rm1, rm2);
    if (qcmnot < 0.0) return {rmnot, 2};

    for (int icount = 1;; ++icount) {
        const double rm = rmmin + rng.uniform() * drm;

        const double pcm = twoBodyP(ww, rm, rm3);
        if (pcm < 0.0) return {rmnot, 3};
        const double qcm = twoBodyP(rm, rm1, rm2);
        if (qcm < 0.0) return {rmnot, 2};

        const double pratio = pcm / pcmnot;
        const double qratio = qcm / qcmnot;
        const double gamma =
            gammanot * (rmnot / rm) * std::pow(qratio, 2.0 * lvalue + 1.0);
        double bwvalue = scale_ * qratio * pratio;
        if (gammanot < 9998.0) {  // otherwise: pure phase space
            const double d1 = rm * rm - rmnot * rmnot;
            bwvalue *= rm * rmnot * gammanot * std::pow(qratio, 2.0 * lvalue) /
                       (d1 * d1 + (rmnot * gamma) * (rmnot * gamma));
        }

        const double ytest = rng.uniform() * bwmaxval_;
        if (bwvalue > bwmaxval_) bwmaxval_ = bwvalue;
        if (ytest <= bwvalue) return {rm, 0};
        if (icount >= 1000) {
            std::printf(" Inefficient BW shape %d %g %g %g %g\n", icount,
                        rmnot, rm1, rm2, rm3);
            return {rm, 0};  // the Fortran also accepts after 1000 tries
        }
    }
}

}  // namespace mcgen
