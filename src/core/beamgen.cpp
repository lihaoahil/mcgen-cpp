#include "core/beamgen.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace mcgen {

Bremsstrahlung::Bremsstrahlung(double pbeammin, double pbeammax)
    : pbeammin_(pbeammin) {
    pendpoint_ = pbeammax / 0.95;  // endpoint specific to the CLAS tagger
    pdelta_ = pbeammax - pbeammin;
    if (pdelta_ <= 0.0)
        throw std::runtime_error("bremsstrahlung: pbeammax <= pbeammin");
    ratio_ = 2.0 * pendpoint_ / 0.000511;
    const double xx = pbeammin / pendpoint_;
    scale_ = (1.0 - xx + 0.75 * xx * xx) *
             (std::log(ratio_ * (1.0 / xx - 1.0)) - 0.5) / pbeammin;
}

double Bremsstrahlung::sample(Rng& rng) const {
    for (;;) {
        const double ppbeam = pbeammin_ + rng.uniform() * pdelta_;
        const double xx = ppbeam / pendpoint_;
        const double factor = (1.0 - xx + 0.75 * xx * xx) *
                              (std::log(ratio_ * (1.0 / xx - 1.0)) - 0.5) /
                              ppbeam;
        if (rng.uniform() <= factor / scale_) return ppbeam;
    }
}

bool BremSculpt::keep(Rng& rng, double ppbeam) {
    // Piecewise-linear GlueX intensity shape (values from the Fortran).
    const double slope1 = 0.15, slope2 = 0.85;
    const double xlo = 3.75, xmid = 7.6, xnot = 8.8, xhi = 11.8;
    double fvalue;
    if (ppbeam <= xmid) {
        fvalue = slope1 * (ppbeam - xlo);
    } else if (ppbeam <= xnot) {
        fvalue = slope1 * (xmid - xlo) + slope2 * (ppbeam - xmid);
    } else if (ppbeam <= xhi) {
        fvalue = 0.27 * (slope1 * (xmid - xlo) + slope2 * (xnot - xmid));
    } else {
        fvalue = 0.0;
    }
    if (fvalue < 0.0) fvalue = 0.0;
    const double ftest = rng.uniform() * fvaluemax_;
    if (fvalue > fvaluemax_) {
        fvaluemax_ = fvalue;
        std::printf(" BREMSCULPT: New largest value of fvalue %g\n", fvaluemax_);
    }
    return ftest <= fvalue;
}

FluxTable::FluxTable(const std::string& path) {
    std::printf(" Reading photon energy distribution from: %s\n", path.c_str());
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open flux file: " + path);
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        double e, c;
        if (!(iss >> e >> c)) continue;
        egamma_.push_back(e - 0.010 / 2.0);  // bin low edge, not center
        counts_.push_back(c);
        if (c > fvaluemax_) fvaluemax_ = c;
    }
    if (egamma_.size() < 2)
        throw std::runtime_error("flux file too short: " + path);

    nmin_ = 0;
    while (nmin_ < counts_.size() && counts_[nmin_] == 0.0) ++nmin_;
    nmax_ = counts_.size() - 1;
    while (nmax_ > 0 && counts_[nmax_] == 0.0) --nmax_;
    std::printf(" First non-0 channel %zu %10.2f %10.2e\n", nmin_ + 1,
                egamma_[nmin_], counts_[nmin_]);
    std::printf(" Last  non-0 channel %zu %10.2f %10.2e\n", nmax_ + 1,
                egamma_[nmax_], counts_[nmax_]);
    std::printf(" Highest flux number %10.2e\n", fvaluemax_);

    const double deltae = std::abs(egamma_[nmin_ + 1] - egamma_[nmin_]);
    if (std::abs(deltae - 0.010) > 0.001)
        throw std::runtime_error(
            "flux file energy bin widths appear not to be 10 MeV");
}

bool FluxTable::keep(Rng& rng, double ppbeam) const {
    if (ppbeam < egamma_[nmin_] || ppbeam > egamma_[nmax_]) return false;
    // Find the first bin edge >= ppbeam and test against the bin below it,
    // exactly as the Fortran does.
    for (size_t n = nmin_; n <= nmax_; ++n) {
        if (egamma_[n] >= ppbeam) {
            if (n == 0) return false;
            return rng.uniform() * fvaluemax_ <= counts_[n - 1];
        }
    }
    return false;
}

Vec3 deuteronMomentum(Rng& rng) {
    constexpr double kAlpha = 0.2316;    // 1/fm
    constexpr double kHbarC = 0.197337;  // GeV fm
    constexpr double kPmaxGen = 0.40;    // GeV/c
    const double alphahc = kAlpha * kHbarC;
    const double betahc = alphahc * 5.98;
    const double a2 = alphahc * alphahc, b2 = betahc * betahc;
    const double phimax = (1.0 / a2 - 1.0 / b2) * (1.0 / a2 - 1.0 / b2);
    for (;;) {
        const double px = (2.0 * rng.uniform() - 1.0) * kPmaxGen;
        const double py = (2.0 * rng.uniform() - 1.0) * kPmaxGen;
        const double pz = (2.0 * rng.uniform() - 1.0) * kPmaxGen;
        const double p2 = px * px + py * py + pz * pz;
        const double phi = (1.0 / (p2 + a2) - 1.0 / (p2 + b2)) *
                           (1.0 / (p2 + a2) - 1.0 / (p2 + b2));
        if (rng.uniform() * phimax <= phi) return {px, py, pz};
    }
}

}  // namespace mcgen
