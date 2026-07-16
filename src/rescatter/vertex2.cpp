#include "rescatter/vertex2.h"

#include <vector>

#include "core/nbody.h"
#include "rescatter/particles.h"

namespace mcgen::rescatter {

bool generateVertex2(Rng& rng, const FourVec& lambdabar, Vertex2Result& out) {
    FourVec target2{0.0, 0.0, 0.0, kMassProton, kMassProton};
    const FourVec p0 = epcm(lambdabar, target2);

    std::vector<double> masses{kMassKaonCh, kMassPionCh, kMassPionCh,
                                kMassPion0};
    std::vector<FourVec> daughters;
    if (!kin45body(rng, p0, masses, daughters)) return false;

    out.kaonPlus = daughters[0];
    out.piPlus = daughters[1];
    out.piMinus = daughters[2];
    out.pi0 = daughters[3];
    return true;
}

}  // namespace mcgen::rescatter
