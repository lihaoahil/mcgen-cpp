#include "rescatter/writer_truth.h"

namespace mcgen::rescatter {

void writeTruthHeader(std::FILE* out) {
    std::fprintf(out,
                 "event,beamE,lambdabar_e,lambdabar_px,lambdabar_py,"
                 "lambdabar_pz,vtx1_x,vtx1_y,vtx1_z,vtx2_x,vtx2_y,vtx2_z\n");
}

void writeTruthRow(std::FILE* out, long eventIndex, const Truth& t) {
    std::fprintf(
        out, "%ld,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
        eventIndex, t.beamEnergy, t.lambdabar.e, t.lambdabar.px,
        t.lambdabar.py, t.lambdabar.pz, t.vertex1[0], t.vertex1[1],
        t.vertex1[2], t.vertex2[0], t.vertex2[1], t.vertex2[2]);
}

}  // namespace mcgen::rescatter
