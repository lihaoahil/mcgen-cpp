#pragma once
// The MC_GEN event generator: main loop, sequential decay chain, and the
// GlueX ASCII event writer.  Ported from the main program of mc_gen.F.

#include <cstdio>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "core/beamgen.h"
#include "core/config.h"
#include "core/dynamics.h"
#include "core/fourvector.h"
#include "core/lineshape.h"
#include "core/rng.h"
#include "core/targetgeom.h"

namespace mcgen {

// One entry of the Fortran PLIST(24,*) table.
struct Particle {
    FourVec p;         // momentum/mass/energy at creation (lab)
    int lund = 0;      // signed Lund ID (999 = initial-state pseudo-particle)
    Vec3 rCreate{};    // creation point (cm)
    Vec3 pol{};        // polarization vector
    FourVec pDecay;    // momentum at decay point (after smearing)
    Vec3 rDecay{};     // decay point (cm)
    double ctau = 0;   // cm
    double charge = 0;
    int index = 0;     // 1-based position in the decay list
    int iparent = 0;   // 1-based parent position (0 for the pseudo-particle)
};

using Event = std::vector<Particle>;  // [0] is the pseudo-particle

// Optional per-event analyzer (mc_analyze); returns true when the event
// passes its cuts ("ikeep").
using Analyzer = std::function<bool(const Event&, double ppbeam)>;

class Generator {
  public:
    Generator(Config cfg, const std::string& defFileName);

    void setAnalyzer(Analyzer a) { analyzer_ = std::move(a); }

    // Run the full generation loop.
    void run();

    Rng& rng() { return rng_; }

  private:
    enum class Status { Done, Abort };  // event chain outcomes

    // One full event attempt; returns the event when the chain completed.
    std::optional<Event> generateEvent();

    double pickBeamMomentum();
    bool decayChain(Event& ev);

    // Decay one particle (plist slot imark); appends products.
    // Returns false to abort the whole event attempt.
    bool decayParticle(Event& ev, size_t imark);

    // Mass selection for the two primary decay products (label 204 block).
    bool pickTwoBodyMasses(int kk, int jj, double& m1, double& m2);

    void writeOutput(Event& ev);
    void openNextOutputFile();
    void printEvent(const Event& ev) const;
    static int lundToGeant(int lund);

    int tableIndexFor(int lund) const;  // -1 if absent; 1-based Fortran index

    Config cfg_;
    Rng rng_;
    Analyzer analyzer_;

    // Beam machinery
    std::optional<Bremsstrahlung> brems_;
    std::optional<FluxTable> flux_;
    BremSculpt sculpt_;

    // Adaptive/stateful dynamics
    ReggePick regge_;
    PhaseSpaceWeight psWeight_;
    PhaseSpaceWeight3 psWeight3_;
    RelativisticLineshape relLineshape_;

    TargetGeom geom_;

    // Event-loop state (Fortran counters and cross-call hacks)
    long ievent_ = 0, numwrite_ = 0, iattempt_ = 0, ikeep_ = 0, icount_ = 0;
    int miterate_ = 0;
    double ppbeam_ = 0.0, wvalue_ = 0.0, ww_ = 0.0;
    double cosinethetacm_ = 0.0, pcmvaluein_ = 0.0;
    double lastPcm_ = 1.0, pcminter_ = 1.0;
    Vec3 pol1_{}, pol2_{};  // persist across decays, as the Fortran locals do
    int imechanism_ = 999;
    bool neutronTarget_ = false;

    // Output stream state
    std::FILE* out_ = nullptr;
    int ifile_ = 0;
    long ieventout_ = 0;
    std::string outputStem_;  // def-file basename without extension
};

}  // namespace mcgen
