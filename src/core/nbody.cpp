#include "core/nbody.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace mcgen {

bool kin3body(Rng& rng, const FourVec& p0, double m1, double m2, double m3,
              FourVec& p1, FourVec& p2, FourVec& p3) {
    const double rm0 = p0.m;
    const double amax23 = rm0 - m1;
    const double amin23 = m2 + m3;
    const double amax12 = rm0 - m3;
    const double amin12 = m1 + m2;

    double rm23 = 0.0;
    for (;;) {
        // Uniform in the squared invariant masses = flat Dalitz density.
        const double rm12 =
            std::sqrt(amin12 * amin12 +
                      (amax12 * amax12 - amin12 * amin12) * rng.uniform());
        rm23 = std::sqrt(amin23 * amin23 +
                         (amax23 * amax23 - amin23 * amin23) * rng.uniform());

        // Local Dalitz limits on m23 for the chosen m12.
        const double e2 = (rm12 * rm12 + m2 * m2 - m1 * m1) / (2.0 * rm12);
        const double e3 = (rm0 * rm0 - rm12 * rm12 - m3 * m3) / (2.0 * rm12);
        const double q2 = std::sqrt(e2 * e2 - m2 * m2);
        const double q3 = std::sqrt(e3 * e3 - m3 * m3);
        const double sum = e2 + e3;
        const double lmax23 = std::sqrt(sum * sum - (q2 - q3) * (q2 - q3));
        const double lmin23 = std::sqrt(sum * sum - (q2 + q3) * (q2 + q3));
        if (rm23 <= lmax23 && rm23 >= lmin23) break;
    }

    // Sequential isotropic decays p0 -> p1 + (23), (23) -> p2 + p3.
    FourVec p23;
    double cth = 0.0, pcm = 0.0;
    if (!kin1r(rng, p0, m1, rm23, p1, p23, cth, pcm)) return false;
    return kin1r(rng, p23, m2, m3, p2, p3, cth, pcm);
}

namespace {

// Shared body of the two quasi-free routines; pickRecoil draws the spectator
// momentum vector (lab frame).
template <typename PickRecoil>
bool quasiFree(Rng& rng, const FourVec& p0, double m1, double m2, double m3,
               PickRecoil&& pickRecoil, FourVec& p1, FourVec& p2, FourVec& p3,
               QuasiFreeResult& qf) {
    if (p0.m - m1 - m2 - m3 < 0.0) {
        std::printf(" No way, no how; masses don't add up: %g %g %g %g\n",
                    p0.m, m1, m2, m3);
        return false;
    }

    FourVec p12;
    for (int imess = 0; imess <= 1000; ++imess) {
        const Vec3 pf = pickRecoil(rng);
        qf.prec = std::sqrt(pf[0] * pf[0] + pf[1] * pf[1] + pf[2] * pf[2]);

        const double erecoil = std::sqrt(qf.prec * qf.prec + m3 * m3);
        p3 = FourVec{pf[0], pf[1], pf[2], m3, 0.0};
        p3.updateEnergy();

        p12.px = p0.px - pf[0];
        p12.py = p0.py - pf[1];
        p12.pz = p0.pz - pf[2];
        const double e12 = p0.etot() - erecoil;
        const double rad = e12 * e12 - p12.p2();
        qf.qfmass = rad > 0.0 ? std::sqrt(rad) : 0.0;
        p12.m = qf.qfmass;
        p12.e = e12;

        if (qf.qfmass >= m1 + m2) {
            double cth = 0.0, pcm = 0.0;
            return kin1r(rng, p12, m1, m2, p1, p2, cth, pcm);
        }
    }
    return false;  // off-shell mass consistently below threshold
}

}  // namespace

bool kin3bodyqf(Rng& rng, const FourVec& p0, double m1, double m2, double m3,
                double pfermi, FourVec& p1, FourVec& p2, FourVec& p3,
                QuasiFreeResult& qf) {
    auto pick = [pfermi](Rng& r) -> Vec3 {
        for (;;) {  // uniform inside the Fermi sphere
            const double px = (2.0 * r.uniform() - 1.0) * pfermi;
            const double py = (2.0 * r.uniform() - 1.0) * pfermi;
            const double pz = (2.0 * r.uniform() - 1.0) * pfermi;
            if (px * px + py * py + pz * pz <= pfermi * pfermi)
                return {px, py, pz};
        }
    };
    return quasiFree(rng, p0, m1, m2, m3, pick, p1, p2, p3, qf);
}

bool kin3bodyhulthen(Rng& rng, const FourVec& p0, double m1, double m2,
                     double m3, FourVec& p1, FourVec& p2, FourVec& p3,
                     QuasiFreeResult& qf) {
    constexpr double kAlpha = 0.2316;    // 1/fm
    constexpr double kHbarC = 0.197337;  // GeV fm
    constexpr double kPmaxGen = 0.40;    // GeV/c
    const double alphahc = kAlpha * kHbarC;
    const double betahc = alphahc * 5.98;
    // The Fortran adaptively tracked the density maximum starting from an
    // overestimate; the Hulthen density peaks at p=0, so use phi(0) exactly.
    const double a2 = alphahc * alphahc, b2 = betahc * betahc;
    const double phimax = (1.0 / a2 - 1.0 / b2) * (1.0 / a2 - 1.0 / b2);

    auto pick = [=](Rng& r) -> Vec3 {
        for (;;) {  // rejection-sample the Hulthen momentum density
            const double px = (2.0 * r.uniform() - 1.0) * kPmaxGen;
            const double py = (2.0 * r.uniform() - 1.0) * kPmaxGen;
            const double pz = (2.0 * r.uniform() - 1.0) * kPmaxGen;
            const double p2 = px * px + py * py + pz * pz;
            const double phi =
                (1.0 / (p2 + a2) - 1.0 / (p2 + b2)) *
                (1.0 / (p2 + a2) - 1.0 / (p2 + b2));
            if (r.uniform() * phimax <= phi) return {px, py, pz};
        }
    };
    return quasiFree(rng, p0, m1, m2, m3, pick, p1, p2, p3, qf);
}

namespace {
// PDK: two-body breakup momentum for parent mass a into masses b, c.
double pdk(double a, double b, double c) {
    const double a2 = a * a, b2 = b * b, c2 = c * c;
    return 0.5 * std::sqrt(a2 + (b2 - c2) * (b2 - c2) / a2 - 2.0 * (b2 + c2));
}

// ROTES2: rotate particle j of pcm by the two GENBOD Euler rotations.
void rotes2(double c, double s, double c2, double s2, FourVec& p) {
    const double sa = p.px;
    const double sb = p.py;
    const double a = sa * c - sb * s;
    p.py = sa * s + sb * c;
    const double b = p.pz;
    p.px = a * c2 - b * s2;
    p.pz = a * s2 + b * c2;
}
}  // namespace

double genbod(Rng& rng, double tecm, const std::vector<double>& masses,
              std::vector<FourVec>& pcm) {
    const int nt = static_cast<int>(masses.size());
    if (nt < 2 || nt > 18) return -1.0;

    double tm = 0.0;
    std::vector<double> sm(static_cast<size_t>(nt));
    for (int i = 0; i < nt; ++i) {
        tm += masses[static_cast<size_t>(i)];
        sm[static_cast<size_t>(i)] = tm;
    }
    const double tecmtm = tecm - tm;
    if (tecmtm <= 0.0) return -1.0;

    // Maximum weight for the constant-cross-section mode (KGENEV=1).
    double wtmax = 1.0;
    {
        double emmax = tecmtm + masses[0];
        double emmin = 0.0;
        for (int i = 1; i < nt; ++i) {
            emmin += masses[static_cast<size_t>(i) - 1];
            emmax += masses[static_cast<size_t>(i)];
            wtmax *= pdk(emmax, emmin, masses[static_cast<size_t>(i)]);
        }
    }
    const double wtmaxq = 1.0 / wtmax;

    // Effective masses: nt-2 ordered uniforms partition the available energy.
    std::vector<double> rno(static_cast<size_t>(3 * nt - 4));
    for (auto& r : rno) r = rng.uniform();
    std::sort(rno.begin(), rno.begin() + (nt - 2));

    std::vector<double> emm(static_cast<size_t>(nt));
    emm[0] = masses[0];
    emm[static_cast<size_t>(nt) - 1] = tecm;
    for (int j = 1; j < nt - 1; ++j)
        emm[static_cast<size_t>(j)] =
            rno[static_cast<size_t>(j) - 1] * tecmtm + sm[static_cast<size_t>(j)];

    double wt = wtmaxq;
    std::vector<double> pd(static_cast<size_t>(nt) - 1);
    for (int i = 0; i < nt - 1; ++i) {
        pd[static_cast<size_t>(i)] = pdk(emm[static_cast<size_t>(i) + 1],
                                         emm[static_cast<size_t>(i)],
                                         masses[static_cast<size_t>(i) + 1]);
        wt *= pd[static_cast<size_t>(i)];
    }

    // Raubold-Lynch construction.
    pcm.assign(static_cast<size_t>(nt), FourVec{});
    for (int i = 0; i < nt; ++i) pcm[static_cast<size_t>(i)].m = masses[static_cast<size_t>(i)];
    pcm[0].py = pd[0];

    int ir = nt - 2;
    for (int i = 1; i < nt; ++i) {
        pcm[static_cast<size_t>(i)].px = 0.0;
        pcm[static_cast<size_t>(i)].py = -pd[static_cast<size_t>(i) - 1];
        pcm[static_cast<size_t>(i)].pz = 0.0;
        const double bang = 2.0 * M_PI * rno[static_cast<size_t>(ir++)];
        const double cb = std::cos(bang), sb = std::sin(bang);
        const double c = 2.0 * rno[static_cast<size_t>(ir++)] - 1.0;
        const double s = std::sqrt(1.0 - c * c);
        if (i < nt - 1) {
            const double esys = std::sqrt(pd[static_cast<size_t>(i)] * pd[static_cast<size_t>(i)] +
                                          emm[static_cast<size_t>(i)] * emm[static_cast<size_t>(i)]);
            const double beta = pd[static_cast<size_t>(i)] / esys;
            const double gama = esys / emm[static_cast<size_t>(i)];
            for (int j = 0; j <= i; ++j) {
                auto& p = pcm[static_cast<size_t>(j)];
                p.updateEnergy();
                rotes2(c, s, cb, sb, p);
                p.py = gama * (p.py + beta * p.e);
            }
        } else {
            for (int j = 0; j <= i; ++j) {
                auto& p = pcm[static_cast<size_t>(j)];
                p.updateEnergy();
                rotes2(c, s, cb, sb, p);
            }
        }
    }
    for (auto& p : pcm) p.updateEnergy();
    return wt;
}

bool kin45body(Rng& rng, const FourVec& p0, const std::vector<double>& masses,
               std::vector<FourVec>& out) {
    double rmtot = 0.0;
    for (double m : masses) rmtot += m;
    if (rmtot >= p0.m) return false;

    std::vector<FourVec> rest;
    if (genbod(rng, p0.m, masses, rest) < 0.0) return false;

    const double energy = p0.etot();
    const Vec3 beta{-p0.px / energy, -p0.py / energy, -p0.pz / energy};
    out.clear();
    out.reserve(rest.size());
    for (const auto& p : rest) out.push_back(boost(p, beta));
    return true;
}

}  // namespace mcgen
