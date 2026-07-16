// mc_gen_rescatter: two-vertex p-Lambdabar rescattering generator.
//   vertex 1: gamma + p1(at rest) -> p + Lambda + Lambdabar
//   vertex 2: Lambdabar + p2(at rest) -> K+ pi+ pi- pi0
// v1: pure phase space at both vertices, flat beam energy, vertex 2
// co-located with vertex 1. See STATUS.md in the pantilamb workspace for
// the full list of deliberate v1 simplifications and the roadmap.
//
// Usage: mc_gen_rescatter <nevents> [beamLoGeV=8.0] [beamHiGeV=9.0]
//                          [outPrefix=rescatter] [seed]
#include <cstdio>
#include <cstdlib>
#include <string>

#include "core/rng.h"
#include "core/targetgeom.h"
#include "rescatter/event.h"
#include "rescatter/writer_ascii.h"
#include "rescatter/writer_truth.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr,
                      "usage: %s <nevents> [beamLoGeV=8.0] [beamHiGeV=9.0] "
                      "[outPrefix=rescatter] [seed]\n",
                      argv[0]);
        return 1;
    }
    const long nevents = std::atol(argv[1]);
    const double beamLo = argc > 2 ? std::atof(argv[2]) : 8.0;
    const double beamHi = argc > 3 ? std::atof(argv[3]) : 9.0;
    const std::string prefix = argc > 4 ? argv[4] : "rescatter";
    mcgen::Rng rng = argc > 5
        ? mcgen::Rng(static_cast<std::uint64_t>(std::atoll(argv[5])))
        : mcgen::Rng();

    // Same LH2 target geometry as validation/gen_mech6_tCutOff.def
    // (cylindrical, 2.0 cm diameter, z in [51, 79] cm, 0.2 cm x/y spread).
    mcgen::TargetGeom geom;
    geom.icode = 2;
    geom.size1 = 2.0;
    geom.size2 = 2.0;
    geom.size3 = 51.0;
    geom.size4 = 79.0;
    geom.jcode = 1;
    geom.xsigma = 0.2;
    geom.ysigma = 0.2;

    std::FILE* asciiOut = std::fopen((prefix + ".ascii").c_str(), "w");
    std::FILE* truthOut = std::fopen((prefix + "_truth.csv").c_str(), "w");
    if (!asciiOut || !truthOut) {
        std::fprintf(stderr, "cannot open output files (prefix=%s)\n",
                     prefix.c_str());
        return 1;
    }
    mcgen::rescatter::writeTruthHeader(truthOut);

    long written = 0, attempts = 0;
    while (written < nevents) {
        ++attempts;
        mcgen::Event ev;
        mcgen::rescatter::Truth truth;
        if (!mcgen::rescatter::generateRescatterEvent(rng, geom, beamLo,
                                                       beamHi, ev, truth))
            continue;
        mcgen::rescatter::writeAsciiEvent(asciiOut, ev, truth.beamEnergy);
        mcgen::rescatter::writeTruthRow(truthOut, written, truth);
        ++written;
    }
    std::fclose(asciiOut);
    std::fclose(truthOut);
    std::printf("wrote %ld events (%ld attempts) to %s.ascii / %s_truth.csv\n",
               written, attempts, prefix.c_str(), prefix.c_str());
    return 0;
}
