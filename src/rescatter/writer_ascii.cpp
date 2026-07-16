#include "rescatter/writer_ascii.h"

#include <vector>

#include "core/particleid.h"

namespace mcgen::rescatter {

void writeAsciiEvent(std::FILE* out, const Event& ev, double beamMomentum,
                     int mechanism) {
    const size_t n = ev.size();
    std::vector<int> geantId(n, -1);
    int icount = 0;
    for (size_t i = 1; i < n; ++i) {
        geantId[i] = lundToGeant(ev[i].lund);
        if (geantId[i] > 0) ++icount;
    }

    std::fprintf(out, "%2d\n", icount);
    std::fprintf(out, "%10.6f%6d\n", beamMomentum, mechanism);

    for (size_t i = 1; i < n; ++i) {
        if (geantId[i] <= 0) continue;
        const auto& q = ev[i];
        std::fprintf(out, "%5zu%5d%5d%10.6f%10.6f%10.6f%10.6f\n", i + 1,
                     q.iparent, geantId[i], q.p.e, q.p.px, q.p.py, q.p.pz);
        std::fprintf(out, "           %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n",
                     q.rCreate[0], q.rCreate[1], q.rCreate[2], q.rDecay[0],
                     q.rDecay[1], q.rDecay[2]);
    }
}

}  // namespace mcgen::rescatter
