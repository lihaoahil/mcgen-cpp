#include "core/rng.h"

#include <chrono>

namespace mcgen {

void Rng::seedFromClock() {
    const auto now = std::chrono::steady_clock::now().time_since_epoch();
    const auto ticks = static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(now).count());
    std::seed_seq seq{static_cast<std::uint32_t>(ticks),
                      static_cast<std::uint32_t>(ticks >> 32),
                      static_cast<std::uint32_t>(std::random_device{}())};
    eng_.seed(seq);
}

}  // namespace mcgen
