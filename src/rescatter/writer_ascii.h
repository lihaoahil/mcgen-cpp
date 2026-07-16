#pragma once
// Writes an mcgen::Event in the standard GlueX generator-input ASCII format
// -- the same line layout as MC_GEN's own Generator::writeOutput, factored
// out so both writers agree on the wire format without duplicating it by
// hand. Verified against validation/check_output.py's parser.

#include <cstdio>

#include "core/generator.h"  // mcgen::Event

namespace mcgen::rescatter {

// mechanism is a free-form integer tag carried through unchanged (MC_GEN
// uses it for its dyncode1-selected production mechanism; we have no such
// concept yet, so it defaults to 0 -- "two-vertex rescattering, v1").
void writeAsciiEvent(std::FILE* out, const Event& ev, double beamMomentum,
                     int mechanism = 0);

}  // namespace mcgen::rescatter
