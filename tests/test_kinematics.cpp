// Physics sanity tests for the ported kinematics:
//  - four-momentum conservation in 2-, 3-, 4-, 5-body decays
//  - invariant masses of products match the table masses
//  - Breit-Wigner samplers stay in range and peak near the pole
// Exits nonzero on any failure.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "core/fourvector.h"
#include "core/lineshape.h"
#include "core/nbody.h"
#include "core/rng.h"

namespace {
int failures = 0;

void check(bool ok, const char* what) {
    if (!ok) {
        std::printf("FAIL: %s\n", what);
        ++failures;
    }
}

void checkClose(double a, double b, double tol, const char* what) {
    if (std::abs(a - b) > tol) {
        std::printf("FAIL: %s (%.9g vs %.9g)\n", what, a, b);
        ++failures;
    }
}

mcgen::FourVec sum(const std::vector<mcgen::FourVec>& v) {
    mcgen::FourVec s;
    for (const auto& p : v) {
        s.px += p.px;
        s.py += p.py;
        s.pz += p.pz;
        s.e += p.e;
    }
    return s;
}
}  // namespace

int main() {
    mcgen::Rng rng(12345);

    // A moving parent: 2.5 GeV mass with 3 GeV/c along a skew direction.
    mcgen::FourVec p0{1.0, -0.5, 3.0, 2.5, 0.0};
    p0.updateEnergy();

    // --- 2-body (kin1r) ---
    for (int i = 0; i < 1000; ++i) {
        mcgen::FourVec p1, p2;
        double cth, pcm;
        check(mcgen::kin1r(rng, p0, 1.115683, 0.493677, p1, p2, cth, pcm),
              "kin1r success");
        auto s = sum({p1, p2});
        checkClose(s.px, p0.px, 1e-9, "kin1r px");
        checkClose(s.py, p0.py, 1e-9, "kin1r py");
        checkClose(s.pz, p0.pz, 1e-9, "kin1r pz");
        checkClose(s.e, p0.etot(), 1e-9, "kin1r E");
        checkClose(p1.e * p1.e - p1.p2(), 1.115683 * 1.115683, 1e-9,
                   "kin1r m1 on-shell");
    }

    // --- 3-body flat phase space ---
    for (int i = 0; i < 1000; ++i) {
        mcgen::FourVec p1, p2, p3;
        check(mcgen::kin3body(rng, p0, 0.938272, 0.139570, 0.139570, p1, p2, p3),
              "kin3body success");
        auto s = sum({p1, p2, p3});
        checkClose(s.px, p0.px, 1e-8, "kin3body px");
        checkClose(s.py, p0.py, 1e-8, "kin3body py");
        checkClose(s.pz, p0.pz, 1e-8, "kin3body pz");
        checkClose(s.e, p0.etot(), 1e-8, "kin3body E");
    }

    // --- 4- and 5-body via GENBOD ---
    for (int nb = 4; nb <= 5; ++nb) {
        std::vector<double> masses(static_cast<size_t>(nb), 0.139570);
        masses[0] = 0.938272;
        for (int i = 0; i < 1000; ++i) {
            std::vector<mcgen::FourVec> out;
            check(mcgen::kin45body(rng, p0, masses, out), "kin45body success");
            auto s = sum(out);
            checkClose(s.px, p0.px, 1e-8, "kin45body px");
            checkClose(s.py, p0.py, 1e-8, "kin45body py");
            checkClose(s.pz, p0.pz, 1e-8, "kin45body pz");
            checkClose(s.e, p0.etot(), 1e-8, "kin45body E");
            for (size_t k = 0; k < out.size(); ++k)
                checkClose(out[k].e * out[k].e - out[k].p2(),
                           masses[k] * masses[k], 1e-8, "kin45body on-shell");
        }
    }

    // --- quasi-free: spectator momentum bounded by pfermi ---
    {
        mcgen::FourVec cm{0.0, 0.0, 3.0, 3.3, 0.0};
        cm.updateEnergy();
        for (int i = 0; i < 500; ++i) {
            mcgen::FourVec p1, p2, p3;
            mcgen::QuasiFreeResult qf;
            check(mcgen::kin3bodyqf(rng, cm, 0.493677, 1.115683, 0.938272,
                                    0.250, p1, p2, p3, qf),
                  "kin3bodyqf success");
            check(qf.prec <= 0.250, "kin3bodyqf recoil inside Fermi sphere");
            auto s = sum({p1, p2, p3});
            checkClose(s.px, cm.px, 1e-8, "qf px");
            checkClose(s.py, cm.py, 1e-8, "qf py");
            checkClose(s.pz, cm.pz, 1e-8, "qf pz");
        }
    }

    // --- bwran: mean ~0, half of the mass within |x|<1 (Cauchy) ---
    {
        int inside = 0;
        const int n = 200000;
        for (int i = 0; i < n; ++i) {
            const double x = mcgen::bwran(rng);
            check(x >= -30.0 && x <= 30.0, "bwran range");
            if (std::abs(x) < 1.0) ++inside;
        }
        // P(|x|<1) for Cauchy truncated to [-30,30] = atan(1)/atan(30) ~ 0.5109
        checkClose(inside / double(n), std::atan(1.0) / std::atan(30.0), 0.01,
                   "bwran CDF at |x|<1");
    }

    // --- relativistic BW: rho -> pi pi, L=1 ---
    {
        mcgen::RelativisticBWSampler bw(0.770, 0.1491, 0.139570, 0.139570, 1);
        double sum1 = 0.0;
        const int n = 100000;
        for (int i = 0; i < n; ++i) {
            const double m = bw.sample(rng);
            check(m >= 2 * 0.139570, "rbw above threshold");
            sum1 += m;
        }
        const double mean = sum1 / n;
        // Long high-side tail pulls the mean above the pole.
        check(mean > 0.770 && mean < 1.2, "rbw mean plausible");
    }

    if (failures == 0) std::printf("all kinematics tests passed\n");
    return failures == 0 ? 0 : 1;
}
