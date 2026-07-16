#include "core/generator.h"

#include <cmath>
#include <cstring>
#include <filesystem>
#include <stdexcept>

#include "core/nbody.h"

namespace mcgen {

namespace {
constexpr double kMassProton = 0.938272081;
constexpr double kMassLambda = 1.115683;
constexpr double kMassDeuteron = 1.875;
}  // namespace

Generator::Generator(Config cfg, const std::string& defFileName)
    : cfg_(std::move(cfg)) {
    // Output stem: basename of the definition file without its extension
    // (mirrors the histnm construction in the Fortran).
    outputStem_ = std::filesystem::path(defFileName).stem().string();

    // RNG seeding.  Deviation from the Fortran GRNDMQ machinery, documented
    // in the README: 0 = fixed seed (reproducible), 1 = continue from
    // randomseed.num when present, 2/3 = seed from the system clock.
    switch (cfg_.irestart) {
        case 0:
            std::printf(" Repeating standard random number sequence\n");
            break;
        case 1: {
            std::FILE* f = std::fopen("randomseed.num", "r");
            unsigned long long s = Rng::kDefaultSeed;
            if (f) {
                if (std::fscanf(f, "%llu", &s) != 1) s = Rng::kDefaultSeed;
                std::fclose(f);
                std::printf(" Continuing saved random number sequence\n");
            }
            rng_.seed(s);
            break;
        }
        case 2:
        case 3:
            std::printf(" Using the system clock for the random number seed\n");
            rng_.seedFromClock();
            break;
        default:
            throw std::runtime_error("no method selected for RNG start");
    }

    geom_.icode = cfg_.targetGeom;
    geom_.size1 = cfg_.size1;
    geom_.size2 = cfg_.size2;
    geom_.size3 = cfg_.size3;
    geom_.size4 = cfg_.size4;
    geom_.jcode = cfg_.vertexCode;
    geom_.xsigma = cfg_.xsigma;
    geom_.ysigma = cfg_.ysigma;

    if (cfg_.iprule == 1 || cfg_.iprule == 2)
        brems_.emplace(cfg_.pbeamMin, cfg_.pbeamMax);
    if (cfg_.iprule == 3) flux_.emplace(cfg_.fluxFile);

    const auto& rm = cfg_.particles;
    neutronTarget_ = rm[1].mass >= 0.939 && rm[1].mass <= 0.940;
    if (neutronTarget_) std::printf(" Neutron target particle\n");
}

int Generator::tableIndexFor(int lund) const {
    // Fortran scans table entries 4..maxnumpar (1-based); return that index.
    for (size_t i = 3; i < cfg_.particles.size(); ++i)
        if (std::abs(lund) == cfg_.particles[i].lundId)
            return static_cast<int>(i) + 1;
    return -1;
}

double Generator::pickBeamMomentum() {
    const double pdelta = cfg_.pbeamMax - cfg_.pbeamMin;
    switch (cfg_.iprule) {
        case 0:
            return pdelta > 0.0 ? cfg_.pbeamMin + rng_.uniform() * pdelta
                                : cfg_.pbeamMin;
        case 1:
            return brems_->sample(rng_);
        case 2:
            for (;;) {
                const double p = brems_->sample(rng_);
                if (sculpt_.keep(rng_, p)) return p;
            }
        case 3:
            for (;;) {
                const double p = pdelta > 0.0
                                     ? cfg_.pbeamMin + rng_.uniform() * pdelta
                                     : cfg_.pbeamMin;
                if (flux_->keep(rng_, p)) return p;
            }
        default:
            throw std::runtime_error(
                "unknown rule requested for selecting beam momentum");
    }
}

bool Generator::pickTwoBodyMasses(int kk, int jj, double& m1, double& m2) {
    // Fortran label 204: choose the rest masses of the two primary decay
    // products, honoring widths / lineshapes.  kk, jj are 1-based table
    // indices.  Returns false when the event must be aborted.
    for (int itwo = 1; itwo <= 2; ++itwo) {
        const int ix = (itwo == 1) ? kk : jj;
        const auto& px = cfg_.particles[static_cast<size_t>(ix) - 1];
        double& mass = (itwo == 1) ? m1 : m2;

        if (px.lundId == 100) {
            // "Box": uniform mass between rm and rm + 2*width.
            const double rmlo = px.mass;
            const double rmhi = px.mass + 2.0 * px.halfWidth;
            mass = rmlo + rng_.uniform() * (rmhi - rmlo);
        } else if (px.halfWidth == 0.0) {
            mass = px.mass;
            if (px.lundId == 60) mass += kMassProton + 0.001;  // magic Lambda hack
        } else if (px.ctau == -1.0) {
            throw std::runtime_error(
                "Lambda(1405) lineshapepick option is not implemented "
                "(also unimplemented in the current Fortran)");
        } else if (px.halfWidth > 0.0) {
            // Pure Lorentzian via bwran; retry on negative mass.
            for (;;) {
                mass = px.mass + px.halfWidth * bwran(rng_);
                if (mass >= 0.0) break;
            }
        } else {
            // Negative width: relativistic BW lineshape with phase space.
            const double gamval = -2.0 * px.halfWidth;
            int lvalue = 0;
            if (px.ctau == -11.0) lvalue = 1;
            else if (px.ctau == -12.0) lvalue = 2;

            // Decay-product masses of the first 2-body mode (GETMASSES).
            if (px.decays.empty() || px.decays[0].nbody != 2)
                throw std::runtime_error(
                    "decay to BW lineshape not defined for this kinematics");
            const int idx1 = tableIndexFor(px.decays[0].products[0]);
            const int idx2 = tableIndexFor(px.decays[0].products[1]);
            if (idx1 < 0 || idx2 < 0)
                throw std::runtime_error(
                    "failed to decode mass values of BW decay products");
            const double rm1 = cfg_.particles[static_cast<size_t>(idx1) - 1].mass;
            const double rm2 = cfg_.particles[static_cast<size_t>(idx2) - 1].mass;

            double rm3 = kMassProton;  // hack, as in the Fortran
            if (px.lundId == 96) rm3 = kMassLambda;

            const auto res = relLineshape_.sample(rng_, px.mass, gamval, rm1,
                                                  rm2, lvalue, rm3, ww_);
            mass = res.mass;
        }
    }
    return true;
}

bool Generator::decayParticle(Event& ev, size_t imark) {
    auto& parent = ev[imark];
    const int lund = parent.lund;

    int iindex;  // 1-based table index
    if (lund == 999) {
        iindex = 3;
    } else {
        iindex = tableIndexFor(lund);
        if (iindex < 0) {
            std::printf(" Could not find particle in table %d\n", lund);
            throw std::runtime_error("unknown particle in decay chain");
        }
    }
    const auto& pdef = cfg_.particles[static_cast<size_t>(iindex) - 1];

    FourVec pv = parent.p;

    // Propagate to the decay point and apply instrumental smearing.
    const double ctau = pdef.ctau > 0.0 ? pdef.ctau : 0.0;
    double zpar = pdef.charge;
    if (lund < 0) zpar = -zpar;
    Vec3 r1 = parent.rCreate;
    FourVec pvsmear = pv;
    if (lund != 999) {
        r1 = propagate(rng_, pv, parent.rCreate, ctau);
        pvsmear = gsmear(rng_, pv, lund, cfg_.dpOverPOut);
    }
    parent.pDecay = pvsmear;
    parent.rDecay = r1;
    parent.ctau = ctau;
    parent.charge = zpar;

    const int ndec = static_cast<int>(pdef.decays.size()) *
                     (pdef.decaysInOrder ? -1 : 1);
    if (pdef.decays.empty()) return true;  // stable

    // Pick a decay mode.
    int mode;  // 1-based
    if (std::abs(ndec) == 1 || ndec < 0) {
        mode = 1;
    } else {
        const double which = rng_.uniform();
        double rsum = 0.0;
        mode = 0;
        for (int m = 1; m <= std::abs(ndec); ++m) {
            const double bran =
                pdef.decays[static_cast<size_t>(m) - 1].fraction;
            if (which >= rsum && which <= rsum + bran) {
                mode = m;
                break;
            }
            rsum += bran;
        }
        if (mode == 0) {
            std::printf(" Failed to select a decay mode\n"
                        " Is the sum of decay branches less than 1.0?\n");
            printEvent(ev);
            return false;  // try another event
        }
    }

    // Label 170: mode retry loop (for decaysInOrder / threshold cases).
    for (;;) {
        const auto& dm = pdef.decays[static_cast<size_t>(mode) - 1];
        const int nbody = dm.nbody;
        int ndynamics = dm.dynamics;
        const double dyncode1 = dm.dyncode1;
        const double dyncode2 = dm.dyncode2;
        const double dyncode3 = dm.dyncode3;

        int ids[5] = {dm.products[0], dm.products[1], dm.products[2],
                      dm.products[3], dm.products[4]};
        if (lund < 0)
            for (int& id : ids) id = -id;

        int tidx[5] = {0, 0, 0, 0, 0};  // 1-based table indices of products
        for (int k = 0; k < nbody; ++k) {
            tidx[k] = tableIndexFor(ids[k]);
            if (tidx[k] < 0) {
                std::printf(" Could not find decay par%d in table %d\n", k + 1,
                            ids[k]);
                throw std::runtime_error("unknown decay product");
            }
        }
        const int kk = tidx[0], jj = tidx[1];

        // Mass selection with the miterate retry logic (labels 204/1306).
        double m[5] = {0, 0, 0, 0, 0};
        bool iwidthflag = false;
        for (;;) {
            if (!pickTwoBodyMasses(kk, jj, m[0], m[1])) return false;
            {
                const auto& pk = cfg_.particles[static_cast<size_t>(kk) - 1];
                const auto& pj = cfg_.particles[static_cast<size_t>(jj) - 1];
                if (pk.halfWidth != 0.0 || pj.halfWidth != 0.0)
                    iwidthflag = true;
            }
            for (int k = 2; k < nbody; ++k) {
                const auto& pd = cfg_.particles[static_cast<size_t>(tidx[k]) - 1];
                if (pd.halfWidth == 0.0) {
                    m[k] = pd.mass;
                } else {
                    m[k] = pd.mass + pd.halfWidth * bwran(rng_);
                    iwidthflag = true;
                }
            }

            double msum = 0.0;
            for (int k = 0; k < nbody; ++k) msum += m[k];
            if (pv.m > msum) break;  // masses OK

            ++miterate_;
            ++iattempt_;
            if (miterate_ == 1) continue;      // re-pick masses once, cheaply
            if (miterate_ >= 10) {             // too many iterations: give up
                std::printf(" Give up on this event %d\n", miterate_);
                miterate_ = 0;
                return false;
            }
            if (iwidthflag) return false;      // retry whole event
            // No widths anywhere: try the next decay mode, if any.
            if (mode < std::abs(ndec)) {
                ++mode;
                goto next_mode;
            }
            return false;
        }

        {
            // Now perform the decay.
            FourVec p1, p2, p3, p4, p5;
            Vec3 pol3{}, pol4{}, pol5{};
            const int id1 = ids[0], id2 = ids[1];

            if (nbody == 2) {
                const Vec3 polv = parent.pol;

                if (ndynamics == 1) {
                    if (pv.p2() == 0.0) {
                        ndynamics = 0;
                        std::printf(" Zero momentum PV decay\n");
                    }
                    const double polnorm2 = polv[0] * polv[0] +
                                            polv[1] * polv[1] +
                                            polv[2] * polv[2];
                    if (polnorm2 == 0.0) ndynamics = 0;
                }

                // Hack: save kaon momentum in the CM frame for later 3-body
                // phase-space weighting.
                pcminter_ = (id1 == +18) ? lastPcm_ : 1.0;

                double cthcm = 0.0;
                bool ok = true;
                switch (ndynamics) {
                    case 0: {
                        double pcm = 0.0;
                        ok = kin1r(rng_, pv, m[0], m[1], p1, p2, cthcm, pcm);
                        if (ok) lastPcm_ = pcm;
                        polarize(pv, p1, p2, id1, id2, pol1_, pol2_, cthcm);
                        break;
                    }
                    case 1:
                        ok = weakDecay(rng_, pv, m[0], m[1], polv, p1, p2,
                                       pol1_, pol2_, cthcm, lund);
                        break;
                    case 2: {
                        double pcm = 0.0;
                        int isum = 0;
                        for (;;) {
                            ok = kin1r(rng_, pv, m[0], m[1], p1, p2, cthcm, pcm);
                            if (p2.p() <= 0.500) break;
                            if (++isum > 500) {
                                std::printf(" Ack! No luck generating small "
                                            "recoil events\n");
                                isum = 0;
                            }
                        }
                        lastPcm_ = pcm;
                        polarize(pv, p1, p2, id1, id2, pol1_, pol2_, cthcm);
                        break;
                    }
                    case 3: {
                        double cth = cosinethetacm_;
                        ok = productStuffit(rng_, pv, m[0], m[1], p1, p2, cth);
                        pol1_ = {0, 0, 0};
                        pol2_ = {0, 0, 0};
                        break;
                    }
                    case 4: {
                        double pcm = 0.0;
                        ok = kin1r(rng_, pv, m[0], m[1], p1, p2, cthcm, pcm);
                        if (ok) lastPcm_ = pcm;
                        pol1_ = polv;  // polarization transfer hack
                        pol2_ = polv;
                        break;
                    }
                    case 5: {
                        double pcm = 0.0;
                        ok = kin1r(rng_, pv, m[0], m[1], p1, p2, cthcm, pcm);
                        if (ok) lastPcm_ = pcm;
                        if (!psWeight_.keep(rng_, pv.m, m[0], m[1], 1.0))
                            return false;
                        break;
                    }
                    case 6: {
                        double cth = cosinethetacm_;
                        ok = productStuffit(rng_, pv, m[0], m[1], p1, p2, cth);
                        pol1_ = {0, 0, 0};
                        pol2_ = {0, 0, 0};
                        if (!psWeight_.keep(rng_, pv.m, m[0], m[1],
                                            pcmvaluein_))
                            return false;
                        break;
                    }
                    case 7: {
                        double pcm = 0.0;
                        ok = kin1r(rng_, pv, m[0], m[1], p1, p2, cthcm, pcm);
                        if (ok) lastPcm_ = pcm;
                        if (!psWeight3_.keep(rng_, pv.m, m[0], m[1],
                                             pcminter_, pcmvaluein_, wvalue_))
                            return false;
                        break;
                    }
                    case 8: {
                        // Regge-type decay: recoil FIRST in the .def list.
                        const double bslope1 = dyncode1;
                        const double tcutoff = dyncode3;
                        const double rmrecoil = m[0];
                        const double rmproduce = m[1];
                        if (!regge_.pick(rng_, ppbeam_, bslope1, tcutoff,
                                         0.001, kMassProton, rmproduce,
                                         rmrecoil, cosinethetacm_)) {
                            std::printf(" Reggepick error\n");
                            return false;
                        }
                        double cth = cosinethetacm_;
                        ok = productStuffit(rng_, pv, m[0], m[1], p1, p2, cth);
                        polarize(pv, p1, p2, id1, id2, pol1_, pol2_,
                                 cosinethetacm_);
                        break;
                    }
                    default:
                        throw std::runtime_error(
                            "unspecified 2-body decay dynamics " +
                            std::to_string(ndynamics));
                }
                if (!ok) std::printf(" KIN1R/decay is messed (nd=%d)\n",
                                     ndynamics);
            } else if (nbody == 3) {
                if (ndynamics == 0) {
                    if (!kin3body(rng_, pv, m[0], m[1], m[2], p1, p2, p3))
                        return false;
                } else if (ndynamics == 1) {
                    QuasiFreeResult qf;
                    if (!kin3bodyqf(rng_, pv, m[0], m[1], m[2], cfg_.pfermi,
                                    p1, p2, p3, qf)) {
                        std::printf(" QF generation failed: discard event\n");
                        return false;
                    }
                } else if (ndynamics == 2) {
                    QuasiFreeResult qf;
                    if (!kin3bodyhulthen(rng_, pv, m[0], m[1], m[2], p1, p2,
                                         p3, qf)) {
                        std::printf(" QF generation failed: discard event\n");
                        return false;
                    }
                } else if (ndynamics >= 3 && ndynamics <= 6) {
                    // Mechanisms 4..7: 3-body phase space in the CM plus a
                    // Regge-style orientation of the (2,3) pair.
                    const double bslope1 = dyncode1;
                    const double bslope2 = dyncode2;
                    const double bslope3 = dyncode3;
                    for (;;) {  // label 240 rethrow loop
                        FourVec pvtemp;
                        pvtemp.m = pv.m;
                        pvtemp.updateEnergy();
                        if (!kin3body(rng_, pvtemp, m[0], m[1], m[2], p1, p2,
                                      p3))
                            return false;
                        FourVec p23;
                        invarMass(p2, p3, p23);
                        const double rm23 = p23.m;

                        double rmproduce, rmrecoil;
                        int ipick = 0;
                        if (ndynamics == 3) {          // Mechanism 4
                            rmproduce = rm23;
                            rmrecoil = m[0];
                            imechanism_ = 4;
                        } else if (ndynamics == 4) {   // Mechanism 5
                            rmproduce = m[0];
                            rmrecoil = rm23;
                            imechanism_ = 5;
                        } else if (ndynamics == 5) {   // Mechanism 6
                            rmproduce = rm23;
                            rmrecoil = m[0];
                            ipick = 1;
                            imechanism_ = 6;
                        } else {                       // Mechanism 7
                            rmproduce = rm23;
                            rmrecoil = m[0];
                            ipick = 2;
                            imechanism_ = 7;
                        }

                        if (ndynamics == 5 || ndynamics == 6) {
                            const double rmassmin = m[1] + m[2];
                            if (!sculptMass(rng_, rm23, rmassmin, bslope2,
                                            bslope3, ipick))
                                continue;  // re-throw the phase space
                        }

                        if (!regge_.pick(rng_, ppbeam_, bslope1, bslope3,
                                         0.001, kMassProton, rmproduce,
                                         rmrecoil, cosinethetacm_)) {
                            std::printf(" Reggepick error\n");
                            return false;
                        }
                        double cth = cosinethetacm_;
                        if (ndynamics == 4) cth = -cth;  // kinematic trick

                        if (!productRotate(rng_, pv, p1, p2, p3, p23, cth))
                            return false;
                        boostToLab(pv, p1, p2, p3, p23);

                        if (ndynamics == 4) {
                            // Double Regge: test the bottom-vertex t.
                            FourVec ptarget;
                            ptarget.m = kMassProton;
                            ptarget.updateEnergy();
                            double tt, ttmin, ttprime;
                            gett(ptarget, p2, tt, ttmin, ttprime);
                            if (!doubleReggePick(rng_, ppbeam_, kMassProton,
                                                 tt, bslope2))
                                continue;  // re-throw the phase space
                        }
                        break;
                    }
                } else {
                    throw std::runtime_error("unspecified 3-body dynamics " +
                                             std::to_string(ndynamics));
                }
            } else if (nbody == 4 || nbody == 5) {
                std::vector<double> masses(m, m + nbody);
                std::vector<FourVec> out;
                if (ndynamics == 0) {
                    if (!kin45body(rng_, pv, masses, out)) return false;
                } else if (ndynamics == 2) {
                    // Demand small spectator (p1) momentum.
                    int isum = 0;
                    for (;;) {
                        if (!kin45body(rng_, pv, masses, out)) return false;
                        if (out[0].p() <= cfg_.pfermi) break;
                        if (++isum > 500) {
                            std::printf(" Ack! No luck generating small "
                                        "recoil events\n");
                            isum = 0;
                        }
                    }
                } else {
                    throw std::runtime_error("unspecified N-body dynamics " +
                                             std::to_string(ndynamics));
                }
                p1 = out[0];
                p2 = out[1];
                p3 = out[2];
                if (nbody >= 4) p4 = out[3];
                if (nbody == 5) p5 = out[4];
            } else {
                throw std::runtime_error(
                    "too many final-state particles requested");
            }

            // Load the products into the event list.  Their creation point
            // is the parent's decay point.  Snapshot the parent fields:
            // push_back may reallocate and invalidate the reference.
            const Vec3 parentDecay = parent.rDecay;
            const int parentIndex = parent.index;
            const FourVec* prods[5] = {&p1, &p2, &p3, &p4, &p5};
            const Vec3* pols[5] = {&pol1_, &pol2_, &pol3, &pol4, &pol5};
            for (int k = 0; k < nbody; ++k) {
                Particle q;
                q.p = *prods[k];
                q.lund = ids[k];
                q.rCreate = parentDecay;
                q.pol = *pols[k];
                q.index = static_cast<int>(ev.size()) + 1;
                q.iparent = parentIndex;
                ev.push_back(q);
            }
            return true;
        }
    next_mode:;
    }
}

bool Generator::decayChain(Event& ev) {
    for (size_t imark = 0; imark < ev.size(); ++imark) {
        if (!decayParticle(ev, imark)) return false;
        if (ev.size() > 60)
            throw std::runtime_error("decay list overflow (maxnumlis)");
    }
    return true;
}

std::optional<Event> Generator::generateEvent() {
    const auto& rm = cfg_.particles;
    imechanism_ = 999;

    ppbeam_ = pickBeamMomentum();
    ww_ = std::sqrt(2.0 * ppbeam_ * kMassProton + kMassProton * kMassProton);
    cosinethetacm_ = 2.0 * rng_.uniform() - 1.0;

    FourVec beam{0.0, 0.0, ppbeam_, rm[0].mass, 0.0};
    beam.updateEnergy();

    FourVec target;
    target.m = rm[1].mass;
    FourVec ptot, ptotrest;
    if (!neutronTarget_) {
        target.updateEnergy();
        ptot = epcm(beam, target);
        ptotrest = ptot;
        wvalue_ = ptot.m;
    } else {
        const Vec3 pn = deuteronMomentum(rng_);
        FourVec moving{pn[0], pn[1], pn[2], rm[1].mass, 0.0};
        moving.updateEnergy();
        ptot = epcm(beam, moving);
        FourVec deut;
        deut.m = kMassDeuteron;
        deut.updateEnergy();
        ptotrest = epcm(beam, deut);
        wvalue_ = ptotrest.m;
    }

    double pcm2 = wvalue_ * wvalue_ - (rm[0].mass * rm[0].mass +
                                       rm[1].mass * rm[1].mass);
    pcm2 = pcm2 * pcm2 -
           (2.0 * rm[0].mass * rm[1].mass) * (2.0 * rm[0].mass * rm[1].mass);
    if (pcm2 <= 0.0) throw std::runtime_error("Yikes! bad initial kinematics");
    pcmvaluein_ = std::sqrt(pcm2) / (2.0 * wvalue_);

    const Vec3 r0 = getVertex(rng_, geom_);

    Event ev;
    ev.reserve(64);  // maxnumlis; also keeps parent references stable
    Particle pseudo;
    pseudo.p = ptot;
    pseudo.lund = 999;
    pseudo.rCreate = r0;
    pseudo.pDecay = ptot;
    pseudo.rDecay = r0;
    pseudo.charge = 999;
    pseudo.index = 1;
    pseudo.iparent = 0;
    ev.push_back(pseudo);

    if (!decayChain(ev)) return std::nullopt;
    return ev;
}

void Generator::run() {
    const auto& rm = cfg_.particles;
    std::printf(" %s\n", cfg_.title.c_str());
    std::printf(" Generating %ld events (or %ld passing events)\n",
                cfg_.nevent, cfg_.ngevent);
    std::printf(" Beam momentum range [%g, %g] GeV/c, rule %d\n",
                cfg_.pbeamMin, cfg_.pbeamMax, cfg_.iprule);
    std::printf(" Incident particle: %s ; target particle: %s\n",
                rm[0].name.c_str(), rm[1].name.c_str());

    for (;;) {
        miterate_ = 0;
        ++numwrite_;
        ++ievent_;
        ++iattempt_;
        if (ievent_ % 5000 == 0)
            std::printf(" Number of events generated so far: %ld %ld\n",
                        ievent_, ikeep_);
        if (ievent_ > cfg_.nevent) break;
        if (ikeep_ > cfg_.ngevent) break;

        // Retry (Fortran label 105) until this attempt yields an event.
        std::optional<Event> ev;
        for (;;) {
            ev = generateEvent();
            if (ev) break;
            if (ievent_ > cfg_.nevent || ikeep_ > cfg_.ngevent) break;
        }
        if (!ev) break;

        if (ievent_ <= cfg_.nprint) {
            if (ievent_ == 1) std::printf(" The first event(s):\n");
            printEvent(*ev);
        }
        ++icount_;
        if (cfg_.doCuts && analyzer_) {
            if (analyzer_(*ev, ppbeam_)) ++ikeep_;
        }
        if (cfg_.doGluex) writeOutput(*ev);
        if (numwrite_ >= cfg_.ngevent) {
            std::printf(" Hit the limit of number of events to write: %ld %ld\n",
                        numwrite_, cfg_.ngevent);
            break;
        }
    }

    std::printf(" Total number of events generated = %ld\n", ievent_ - 1);
    if (out_) std::fclose(out_);

    // Save a continuation seed (deviation: single-seed scheme).
    std::FILE* f = std::fopen("randomseed.num", "w");
    if (f) {
        std::fprintf(f, "%llu\n",
                     static_cast<unsigned long long>(
                         rng_.uniform() * 9.0e18));
        std::fclose(f);
    }
    std::printf(" Successful end of event generation.\n");
}

int Generator::lundToGeant(int lund) {
    switch (lund) {
        case 1: case -1: return 1;    // photon
        case -7: return 2;            // positron
        case 7: return 3;             // electron
        case 2: return 4;             // neutrino
        case -9: return 5;            // mu+
        case 9: return 6;             // mu-
        case 23: case -23: return 7;  // pi0
        case 17: return 8;            // pi+
        case -17: return 9;           // pi-
        case 18: return 11;           // K+
        case -18: return 12;          // K-
        case 42: return 13;           // neutron
        case 41: return 14;           // proton
        case -41: return 15;          // anti-proton
        case 19: case -19: return 16; // K short
        case 24: return 17;           // eta
        case 57: return 18;           // Lambda
        case -57: return 26;          // anti-Lambda
        case 43: return 19;           // Sigma+
        case 44: case -44: return 20; // Sigma0
        case 45: return 21;           // Sigma-
        case 46: return 23;           // Xi-
        case -46: return 31;          // anti-Xi+
        case 47: return 22;           // Xi0
        case -47: return 30;          // anti-Xi0
        case 33: return 57;           // rho0
        case 27: return 58;           // rho+
        case -27: return 59;          // rho-
        case 34: return 60;           // omega
        case 36: return 61;           // eta prime
        case 35: return 62;           // phi
        case 20: case 21: return 0;   // K*0, K*+ (decayed here, not in GEANT)
        case 91: return 991;
        case 92: return 992;
        case 93: return 993;
        case 94: return 994;
        case 95: return 995;
        case 96: return 996;
        case 97: return 997;
        case 98: return 998;
        case 99: return 999;
        case 100: return 9100;
        case 101: case -101: return 9101;
        default: return -1;
    }
}

void Generator::openNextOutputFile() {
    if (out_) std::fclose(out_);
    ++ifile_;
    char nam[8];
    std::snprintf(nam, sizeof nam, "%04d", ifile_);
    const std::string fname =
        cfg_.prefix + outputStem_ + "_" + nam + ".ascii";
    std::printf(" New Monte Carlo output %s\n", fname.c_str());
    out_ = std::fopen(fname.c_str(), "w");
    if (!out_)
        throw std::runtime_error("cannot open output file: " + fname);
}

void Generator::writeOutput(Event& ev) {
    if (!out_) {
        ieventout_ = 0;
        openNextOutputFile();
    }
    ++ieventout_;
    if (ieventout_ > cfg_.maxOutputPerFile) {
        ieventout_ = 1;
        openNextOutputFile();
    }

    const size_t n = ev.size();
    std::vector<int> geantId(n, -1);
    int icount = 0, iproton = 0;
    for (size_t i = 1; i < n; ++i) {
        geantId[i] = lundToGeant(ev[i].lund);
        if (geantId[i] == -1)
            std::printf(" Unknown Lund-->HDGEANT mapping for Lund ID: %d\n",
                        ev[i].lund);
        if (geantId[i] > 0) ++icount;
        if (geantId[i] == 14) ++iproton;
    }

    // Lambda-anti-Lambda hack: HDDM wants the non-Lambda-decay proton first.
    // "Mechanism 2" (kaon exchange, ID 96 in slot 3) needs a rotation of
    // slots 2..6 (Fortran indices) = ev[1..5].
    if (iproton == 2 && icount == 5 && n > 5 && ev[2].lund == 96) {
        Particle tmp = ev[5];
        int tmpId = geantId[5];
        for (size_t j = 5; j > 1; --j) {
            ev[j] = ev[j - 1];
            geantId[j] = geantId[j - 1];
        }
        ev[1] = tmp;
        geantId[1] = tmpId;
    }

    // Sanity checks: NaN or absurd vertex values abort the event.
    const double toobig = 999.0;  // centimeters
    for (size_t i = 1; i < n; ++i) {
        if (geantId[i] <= 0) continue;
        const auto& q = ev[i];
        const double vals[] = {q.rCreate[0], q.rCreate[1], q.rCreate[2],
                               q.rDecay[0],  q.rDecay[1],  q.rDecay[2],
                               q.p.e,        q.p.px,       q.p.py,
                               q.p.pz};
        for (double v : vals) {
            if (std::isnan(v)) {
                std::printf(" Aborting event: NaN status\n");
                --ieventout_;
                return;
            }
        }
        if (std::abs(q.rCreate[0]) > toobig ||
            std::abs(q.rCreate[1]) > toobig) {
            std::printf(" Aborting event: vertex too large\n");
            --ieventout_;
            return;
        }
    }

    // Event header: particle count, then beam momentum + mechanism.
    std::fprintf(out_, "%2d\n", icount);
    std::fprintf(out_, "%10.6f%6d\n", ppbeam_, imechanism_);

    for (size_t i = 1; i < n; ++i) {
        if (geantId[i] <= 0) continue;
        const auto& q = ev[i];
        double vxd = q.rDecay[0], vyd = q.rDecay[1], vzd = q.rDecay[2];
        if (std::abs(vxd) > toobig) vxd = toobig;
        if (std::abs(vyd) > toobig) vyd = toobig;
        if (std::abs(vzd) > toobig) vzd = toobig;
        std::fprintf(out_, "%5zu%5d%5d%10.6f%10.6f%10.6f%10.6f\n", i + 1,
                     q.iparent, geantId[i], q.p.e, q.p.px, q.p.py, q.p.pz);
        std::fprintf(out_, "           %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n",
                     q.rCreate[0], q.rCreate[1], q.rCreate[2], vxd, vyd, vzd);
    }
}

void Generator::printEvent(const Event& ev) const {
    std::printf("\n Event %7ld; Beam momentum%8.2f\n", ievent_, ppbeam_);
    for (const auto& q : ev) {
        const double beta = q.p.p() / q.p.e;
        const double tof = 4.0 / (beta * 0.3);
        std::printf(" Indicies:  Particle:%3d; Parent:%3d; LUND ID:%4d\n",
                    q.index, q.iparent, q.lund);
        std::printf("       px,py,pz,m,E %10.3f%10.3f%10.3f%10.3f%10.3f"
                    "            beta, tof %10.3f%10.3f\n",
                    q.p.px, q.p.py, q.p.pz, q.p.m, q.p.e, beta, tof);
        std::printf("       X,Y,Z         %10.1e%10.1e%10.1e  Pol %8.4f%8.4f%8.4f\n",
                    q.rCreate[0], q.rCreate[1], q.rCreate[2], q.pol[0],
                    q.pol[1], q.pol[2]);
    }
}

}  // namespace mcgen
