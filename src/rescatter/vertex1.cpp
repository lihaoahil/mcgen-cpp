#include "rescatter/vertex1.h"

#include "core/nbody.h"
#include "rescatter/particles.h"

namespace mcgen::rescatter {

bool generateVertex1(Rng& rng, double beamEnergy, Vertex1Result& out) {
    FourVec beam{0.0, 0.0, beamEnergy, 0.0, beamEnergy};
    FourVec target{0.0, 0.0, 0.0, kMassProton, kMassProton};
    const FourVec p0 = epcm(beam, target);
    return kin3body(rng, p0, kMassProton, kMassLambda, kMassLambda,
                     out.proton, out.lambda, out.lambdabar);
}

}  // namespace mcgen::rescatter
