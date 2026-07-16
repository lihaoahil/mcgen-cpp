#pragma once
// Beam-photon energy selection, ported from mc_gen.F: BREMSSTRAHLUNG,
// BREMSCULPT, BREMSELECT (external GlueX flux file), GETDEUTERONMOMENTUM.

#include <string>
#include <vector>

#include "core/fourvector.h"
#include "core/rng.h"

namespace mcgen {

// IPRULE=1: photon energy from a bremsstrahlung spectrum shape (box approx;
// endpoint = pbeammax/0.95, CLAS-tagger specific, as in the Fortran).
class Bremsstrahlung {
  public:
    Bremsstrahlung(double pbeammin, double pbeammax);
    double sample(Rng& rng) const;

  private:
    double pbeammin_, pdelta_, pendpoint_, ratio_, scale_;
};

// IPRULE=2 add-on: ad hoc piecewise-linear sculpting of the accepted photon
// spectrum (GlueX Lambda-anti-Lambda shape hack from the Fortran).
class BremSculpt {
  public:
    bool keep(Rng& rng, double ppbeam);

  private:
    double fvaluemax_ = 0.0;  // adaptive envelope
};

// IPRULE=3: accept/reject against a tabulated GlueX beam flux file
// (two columns: bin-center energy in GeV, counts; 10 MeV bins).
class FluxTable {
  public:
    explicit FluxTable(const std::string& path);
    bool keep(Rng& rng, double ppbeam) const;

  private:
    std::vector<double> egamma_;  // bin LOW edges after the -5 MeV shift
    std::vector<double> counts_;
    size_t nmin_ = 0, nmax_ = 0;  // first/last non-empty bins
    double fvaluemax_ = 0.0;
};

// GETDEUTERONMOMENTUM: initial-state nucleon momentum from the deuteron
// Hulthen distribution (shared shape with kin3bodyhulthen).
Vec3 deuteronMomentum(Rng& rng);

}  // namespace mcgen
