#pragma once
// Random number facade replacing the Fortran RAN / RNDNOR / RANNOR suite
// (rndnor.F).  The port is statistically equivalent, not bit-compatible:
// std::mt19937_64 replaces the legacy congruential generator, and
// std::normal_distribution replaces the Box-Muller RANNOR.
//
// The Fortran seeded either from a fixed table (irestart=0, reproducible) or
// from the system clock; Rng mirrors that with a fixed default seed and an
// optional clock seed.

#include <cstdint>
#include <random>

namespace mcgen {

class Rng {
  public:
    explicit Rng(std::uint64_t seed = kDefaultSeed) : eng_(seed) {}

    void seed(std::uint64_t s) { eng_.seed(s); }
    void seedFromClock();

    // RAN(ISEED): uniform in [0,1).
    double uniform() { return uni_(eng_); }

    // Uniform in [a,b).
    double uniform(double a, double b) { return a + (b - a) * uni_(eng_); }

    // RANNOR/RNDNOR: standard normal deviate.
    double normal() { return norm_(eng_); }

    // Gaussian with given mean and sigma.
    double normal(double mean, double sigma) { return mean + sigma * norm_(eng_); }

    static constexpr std::uint64_t kDefaultSeed = 20260702ULL;

  private:
    std::mt19937_64 eng_;
    std::uniform_real_distribution<double> uni_{0.0, 1.0};
    std::normal_distribution<double> norm_{0.0, 1.0};
};

}  // namespace mcgen
