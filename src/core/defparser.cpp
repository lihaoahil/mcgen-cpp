// Parser for MC_GEN .def definition files, ported from the input section of
// mc_gen.F (labels 3..40).
//
// Format notes, verified empirically against gfortran's behavior:
//  - Everything before a line starting with "=====" is a comment block; the
//    line right after it is the run title.
//  - Header values are read one per line with Fortran err-retry semantics:
//    a line that fails to parse is skipped and the read retried on the next
//    line.  Trailing "!comments" after values are ignored.
//  - After the header, exactly 13 comment lines are skipped (fixed count in
//    the Fortran).
//  - Particle table entries:  index 'name' mass width zpart ctau lundid
//    ndecays, terminated by index -1.  The width column is (Gamma in MeV);
//    it is stored as Gamma/2 in GeV, sign preserved (negative selects a
//    relativistic Breit-Wigner lineshape).
//  - Each decay mode is a dynamics line (nbody ndyn dyn1 dyn2 dyn3) followed
//    by a product line (nbody Lund IDs + branching fraction).
//
// Intentional deviation from the Fortran: dynamics lines with only 4 values
// get dyncode3 = 0 and a warning.  The Fortran's list-directed read would
// instead spill into the product line and scramble the table (the stale
// gen_lamlambar_mech1.def template trips this), which is never what the
// file author intended.

#include <charconv>
#include <cstdio>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>

#include "core/config.h"

namespace mcgen {

namespace {

// Split a line into whitespace-separated tokens, honoring the fact that a
// token starting with '!' begins a comment in practice (Fortran never reads
// that far because reads stop once their item list is satisfied).
std::vector<std::string> tokenize(const std::string& line) {
    std::vector<std::string> out;
    const size_t n = line.size();
    size_t i = 0;
    while (i < n) {
        while (i < n && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
        if (i >= n) break;
        if (line[i] == '!') break;  // comment to end of line
        std::string tok;
        if (line[i] == '\'') {  // Fortran character constant, may hold spaces
            const size_t close = line.find('\'', i + 1);
            const size_t end = (close == std::string::npos) ? n : close + 1;
            tok = line.substr(i, end - i);
            i = end;
        } else {
            const size_t start = i;
            while (i < n && !std::isspace(static_cast<unsigned char>(line[i])))
                ++i;
            tok = line.substr(start, i - start);
        }
        out.push_back(std::move(tok));
    }
    return out;
}

std::optional<double> toDouble(const std::string& s) {
    try {
        size_t pos = 0;
        // Fortran-style exponents like 1.E32 are accepted by strtod.
        double v = std::stod(s, &pos);
        if (pos != s.size()) return std::nullopt;
        return v;
    } catch (...) {
        return std::nullopt;
    }
}

std::optional<long> toLong(const std::string& s) {
    auto d = toDouble(s);
    if (!d) return std::nullopt;
    return static_cast<long>(*d);
}

// Strip the quotes from a Fortran character constant ('PHOTON' -> PHOTON).
std::string unquote(const std::string& s) {
    if (s.size() >= 2 && s.front() == '\'' && s.back() == '\'')
        return s.substr(1, s.size() - 2);
    return s;
}

class LineReader {
  public:
    explicit LineReader(const std::string& path) : in_(path), path_(path) {
        if (!in_) throw std::runtime_error("cannot open definition file: " + path);
    }

    // Raw next line; throws at EOF (the Fortran would abort similarly).
    std::string nextLine() {
        std::string line;
        if (!std::getline(in_, line))
            throw std::runtime_error(path_ + ": unexpected end of file at line " +
                                     std::to_string(lineno_));
        ++lineno_;
        return line;
    }

    // Fortran err-retry read: skip lines until one yields at least `need`
    // parseable leading numeric tokens; returns the tokens of that line.
    std::vector<std::string> nextNumericLine(size_t need) {
        for (;;) {
            auto toks = tokenize(nextLine());
            if (toks.size() < need) continue;
            bool ok = true;
            for (size_t k = 0; k < need; ++k)
                if (!toDouble(toks[k])) { ok = false; break; }
            if (ok) return toks;
        }
    }

    // read(1,'(a)') equivalent for filename-like fields: first token of the
    // next non-empty line.
    std::string nextStringField() {
        for (;;) {
            auto toks = tokenize(nextLine());
            if (!toks.empty()) return toks[0];
        }
    }

    int lineno() const { return lineno_; }

  private:
    std::ifstream in_;
    std::string path_;
    int lineno_ = 0;
};

}  // namespace

int Config::findByLund(int lundId) const {
    for (size_t i = 0; i < particles.size(); ++i)
        if (particles[i].lundId == std::abs(lundId)) return static_cast<int>(i);
    return -1;
}

Config parseDefFile(const std::string& path) {
    LineReader rd(path);
    Config cfg;

    // Skip the free comment block.
    for (;;) {
        if (rd.nextLine().rfind("=====", 0) == 0) break;
    }
    cfg.title = rd.nextLine();

    auto num = [&rd](size_t n = 1) { return rd.nextNumericLine(n); };

    cfg.nevent = *toLong(num()[0]);
    cfg.ngevent = *toLong(num()[0]);
    cfg.pbeamMin = *toDouble(num()[0]);
    cfg.pbeamMax = *toDouble(num()[0]);
    cfg.iprule = static_cast<int>(*toLong(num()[0]));
    cfg.fluxFile = rd.nextStringField();
    cfg.nprint = static_cast<int>(*toLong(num()[0]));
    cfg.hackFlag = *toLong(num()[0]) == 1;
    cfg.doGluex = *toLong(num()[0]) == 1;
    cfg.doCuts = *toLong(num()[0]) == 1;
    cfg.prefix = rd.nextStringField();
    cfg.maxOutputPerFile = *toLong(num()[0]);
    cfg.dpOverPOut = *toDouble(num()[0]);
    cfg.pfermi = *toDouble(num()[0]);
    cfg.targetGeom = static_cast<int>(*toLong(num()[0]));
    {
        auto t = num(4);
        cfg.size1 = *toDouble(t[0]);
        cfg.size2 = *toDouble(t[1]);
        cfg.size3 = *toDouble(t[2]);
        cfg.size4 = *toDouble(t[3]);
    }
    {
        auto t = num(3);
        cfg.vertexCode = static_cast<int>(*toLong(t[0]));
        cfg.xsigma = *toDouble(t[1]);
        cfg.ysigma = *toDouble(t[2]);
    }
    cfg.irestart = static_cast<int>(*toLong(num()[0]));
    cfg.reserved = *toDouble(num()[0]);

    // Fixed block of 13 comment lines (see mc_gen.F).
    for (int i = 0; i < 13; ++i) rd.nextLine();

    // Particle table.
    for (int i = 1;; ++i) {
        // Particle line: index 'name' mass width zpart ctau lundid ndecays.
        // The name is non-numeric, so nextNumericLine(1) + a shape check.
        std::vector<std::string> toks;
        for (;;) {
            toks = tokenize(rd.nextLine());
            if (toks.size() >= 1 && toDouble(toks[0])) {
                if (*toLong(toks[0]) == -1) break;       // terminator
                if (toks.size() >= 8 && !toDouble(toks[1]) && toDouble(toks[2]))
                    break;
            }
        }
        if (*toLong(toks[0]) == -1) break;

        const long idummy = *toLong(toks[0]);
        if (idummy != i)
            throw std::runtime_error(
                path + ":" + std::to_string(rd.lineno()) +
                ": definition table index mismatch: expected particle " +
                std::to_string(i) + ", file says " + std::to_string(idummy));

        ParticleDef pd;
        pd.name = unquote(toks[1]);
        pd.mass = *toDouble(toks[2]);
        pd.halfWidth = *toDouble(toks[3]) / 2.0 / 1000.0;  // Gamma/2 in GeV
        pd.charge = *toDouble(toks[4]);
        pd.ctau = *toDouble(toks[5]);
        pd.lundId = static_cast<int>(*toLong(toks[6]));
        const long ndecays = *toLong(toks[7]);
        pd.decaysInOrder = ndecays < 0;

        // Beam and target (i <= 2) cannot decay.
        const long nmodes = (i <= 2) ? 0 : std::labs(ndecays);
        for (long j = 0; j < nmodes; ++j) {
            DecayMode dm;
            auto dyn = rd.nextNumericLine(4);
            dm.nbody = static_cast<int>(*toLong(dyn[0]));
            dm.dynamics = static_cast<int>(*toLong(dyn[1]));
            dm.dyncode1 = *toDouble(dyn[2]);
            dm.dyncode2 = *toDouble(dyn[3]);
            if (dyn.size() >= 5 && toDouble(dyn[4])) {
                dm.dyncode3 = *toDouble(dyn[4]);
            } else {
                dm.dyncode3 = 0.0;
                std::fprintf(stderr,
                             "%s:%d: warning: dynamics line for '%s' decay %ld "
                             "has 4 values; dyncode3 set to 0 (the Fortran "
                             "reader expects 5 since the 04-2020 format)\n",
                             path.c_str(), rd.lineno(), pd.name.c_str(), j + 1);
            }
            if (dm.nbody < 2 || dm.nbody > 5)
                throw std::runtime_error(
                    path + ":" + std::to_string(rd.lineno()) + ": N=" +
                    std::to_string(dm.nbody) + " body decay not implemented");

            auto prod = rd.nextNumericLine(static_cast<size_t>(dm.nbody) + 1);
            for (int k = 0; k < dm.nbody; ++k)
                dm.products[static_cast<size_t>(k)] =
                    static_cast<int>(*toLong(prod[static_cast<size_t>(k)]));
            dm.fraction = *toDouble(prod[static_cast<size_t>(dm.nbody)]);
            pd.decays.push_back(dm);
        }
        cfg.particles.push_back(std::move(pd));
    }

    if (cfg.particles.size() < 3)
        throw std::runtime_error(path +
                                 ": particle table needs at least beam, "
                                 "target, and the beam/target pseudo-particle");
    return cfg;
}

}  // namespace mcgen
