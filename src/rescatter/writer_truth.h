#pragma once
// Plain-CSV MC-truth sidecar: the true Lambdabar four-vector and true
// annihilation vertex the ASCII stream can't carry (the Lambdabar is
// consumed, never a final-state particle). Consumed downstream by
// pantilamb/validation scripts to score the -m_p pseudo-particle
// reconstruction strategy against truth.

#include <cstdio>

#include "rescatter/event.h"

namespace mcgen::rescatter {

void writeTruthHeader(std::FILE* out);
void writeTruthRow(std::FILE* out, long eventIndex, const Truth& t);

}  // namespace mcgen::rescatter
