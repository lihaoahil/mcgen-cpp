// MC_GEN C++ port - entry point.
// Usage mirrors the Fortran: mc_gen name_of_definition_file.def
#include <cstdio>
#include <exception>
#include <filesystem>
#include <memory>

#include "core/config.h"
#include "core/generator.h"

#ifdef MCGEN_WITH_ROOT
#include "analysis/rootanalyzer.h"
#endif

int main(int argc, char** argv) {
    std::printf("****************************************************\n");
    std::printf("Carnegie Mellon Monte Carlo Event Generator\n");
    std::printf("  MC_GEN C++20 port (from Fortran v12.92)\n");
    std::printf("****************************************************\n");
    if (argc < 2) {
        std::printf(" Usage: %s name_of_definition_file.def\n", argv[0]);
        return 1;
    }
    try {
        mcgen::Config cfg = mcgen::parseDefFile(argv[1]);
        const bool doCuts = cfg.doCuts;
        mcgen::Generator gen(std::move(cfg), argv[1]);

#ifdef MCGEN_WITH_ROOT
        std::unique_ptr<mcgen::RootAnalyzer> analyzer;
        if (doCuts) {
            // The Fortran wrote <defname>.hbk; the ROOT port writes .root.
            const std::string histFile =
                std::filesystem::path(argv[1]).stem().string() + ".root";
            analyzer = std::make_unique<mcgen::RootAnalyzer>(histFile);
            gen.setAnalyzer(analyzer->asCallback());
        }
#else
        if (doCuts)
            std::printf(" Analyzer requested but built without ROOT; "
                        "events will not be histogrammed.\n");
#endif

        gen.run();

#ifdef MCGEN_WITH_ROOT
        if (analyzer) analyzer->finish();
#endif
    } catch (const std::exception& e) {
        std::fprintf(stderr, "mc_gen: fatal: %s\n", e.what());
        return 1;
    }
    return 0;
}
