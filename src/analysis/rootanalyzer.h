#pragma once
// ROOT-backed port of mc_analyze.F (R. A. Schumacher, CMU): per-event
// analysis and histogramming for the GlueX p-pbar-p, Lambda-anti-Lambda,
// and Ks-Ks-p studies.
//
// HBOOK -> ROOT mapping: histogram IDs become "h<ID>" TH1D/TH2D objects,
// booked identically in three TDirectories (RAWEVENT, ACCEPTED, WEIGHTED)
// of one .root file, mirroring the PAW directory layout.  PAW title markup
// is converted to ROOT TLatex best-effort.
//
// This header keeps ROOT types out of the core: the implementation is
// hidden behind a pimpl.

#include <memory>
#include <string>

#include "core/generator.h"

namespace mcgen {

class RootAnalyzer {
  public:
    // histFile: output .root path (the Fortran wrote <defname>.hbk).
    explicit RootAnalyzer(const std::string& histFile);
    ~RootAnalyzer();

    // Per-event entry point (port of subroutine mc_analyze).  Returns true
    // when the event passes the trigger selection ("ikeep" events).
    bool analyze(const Event& ev, double ppbeam);

    // Port of histogram_finish + analysis_finish: writes the .root file and
    // prints the Lambda anti-Lambda polarization / correlation summary.
    void finish();

    // Adapter for Generator::setAnalyzer.
    Analyzer asCallback() {
        return [this](const Event& ev, double ppbeam) {
            return analyze(ev, ppbeam);
        };
    }

  private:
    struct Impl;

    // The Lambda anti-Lambda (5-track) and Ks Ks p analysis blocks; return
    // false when the event fails an in-block analysis cut (Fortran goto 1000).
    bool lambdaBlock(const Event& ev, double ppbeam, double weight, int iloop,
                     bool iprocess);
    bool ksksBlock(const Event& ev, double ppbeam, double weight, int iloop,
                   bool iprocess, const std::vector<double>& acceptlist);

    std::unique_ptr<Impl> impl_;
};

}  // namespace mcgen
