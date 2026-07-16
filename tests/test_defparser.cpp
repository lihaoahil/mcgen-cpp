// Parses a .def file and dumps the table; exits nonzero on parse failure.
#include <cstdio>

#include "core/config.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s file.def\n", argv[0]);
        return 2;
    }
    try {
        const mcgen::Config cfg = mcgen::parseDefFile(argv[1]);
        std::printf("title: %s\n", cfg.title.c_str());
        std::printf("nevent=%ld ngevent=%ld pbeam=[%g,%g] iprule=%d\n",
                    cfg.nevent, cfg.ngevent, cfg.pbeamMin, cfg.pbeamMax,
                    cfg.iprule);
        std::printf("flux='%s' prefix='%s' perfile=%ld\n", cfg.fluxFile.c_str(),
                    cfg.prefix.c_str(), cfg.maxOutputPerFile);
        std::printf("dpoverp=%g pfermi=%g geom=%d sizes=%g,%g,%g,%g "
                    "vtx=%d(%g,%g) irestart=%d\n",
                    cfg.dpOverPOut, cfg.pfermi, cfg.targetGeom, cfg.size1,
                    cfg.size2, cfg.size3, cfg.size4, cfg.vertexCode,
                    cfg.xsigma, cfg.ysigma, cfg.irestart);
        std::printf("particles: %zu\n", cfg.particles.size());
        for (size_t i = 0; i < cfg.particles.size(); ++i) {
            const auto& p = cfg.particles[i];
            std::printf("%3zu %-10s m=%8.5f G/2=%9.5f q=%+g ctau=%g lund=%d "
                        "modes=%zu%s\n",
                        i + 1, p.name.c_str(), p.mass, p.halfWidth, p.charge,
                        p.ctau, p.lundId, p.decays.size(),
                        p.decaysInOrder ? " [in-order]" : "");
            for (const auto& d : p.decays) {
                std::printf("      nb=%d dyn=%d codes=(%g,%g,%g) -> ", d.nbody,
                            d.dynamics, d.dyncode1, d.dyncode2, d.dyncode3);
                for (int k = 0; k < d.nbody; ++k)
                    std::printf("%d ", d.products[static_cast<size_t>(k)]);
                std::printf(" frac=%g\n", d.fraction);
            }
        }
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "parse error: %s\n", e.what());
        return 1;
    }
}
