#pragma once
// Particle-ID mapping shared by every ASCII event writer (the .def-driven
// Generator and any other event source, e.g. src/rescatter).

namespace mcgen {

// Map an MC_GEN "lund"-table particle ID to the GEANT ID written to the
// GlueX ASCII generator-input stream. Returns 0 for IDs that are
// intentionally not written (decayed internally, e.g. K*), -1 for IDs with
// no known mapping.
int lundToGeant(int lund);

}  // namespace mcgen
