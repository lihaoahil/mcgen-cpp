#pragma once
// PDG masses (GeV) and MC_GEN-table "lund" IDs for the p-Lambdabar
// rescattering final state. See core/particleid.h for the map from these
// lund IDs to the GEANT IDs actually written to the ASCII stream.

namespace mcgen::rescatter {

inline constexpr double kMassProton = 0.938272;
inline constexpr double kMassLambda = 1.115683;  // same for Lambda and Lambdabar
inline constexpr double kMassKaonCh = 0.493677;
inline constexpr double kMassPionCh = 0.139570;
inline constexpr double kMassPion0 = 0.134977;

inline constexpr int kLundProton = 41;
inline constexpr int kLundLambda = 57;
inline constexpr int kLundLambdaBar = -57;
inline constexpr int kLundKaonPlus = 18;
inline constexpr int kLundPiPlus = 17;
inline constexpr int kLundPiMinus = -17;
inline constexpr int kLundPi0 = 23;

}  // namespace mcgen::rescatter
