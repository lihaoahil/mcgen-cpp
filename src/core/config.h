#pragma once
// Configuration ("definition file") data model for MC_GEN, mirroring the
// .def format read by mc_gen.F.  See defparser.cpp for format notes.

#include <array>
#include <string>
#include <vector>

namespace mcgen {

struct DecayMode {
    int nbody = 0;                     // 2..5 bodies
    int dynamics = 0;                  // IDECAY(i,j,7): decay-dynamics code
    double dyncode1 = 0.0;             // BRANCH(i,j,2)
    double dyncode2 = 0.0;             // BRANCH(i,j,3)
    double dyncode3 = 0.0;             // BRANCH(i,j,4)
    std::array<int, 5> products{};     // Lund IDs (negative = antiparticle)
    double fraction = 0.0;             // BRANCH(i,j,1) branching fraction
};

struct ParticleDef {
    std::string name;
    double mass = 0.0;                 // GeV
    double halfWidth = 0.0;            // Gamma/2 in GeV (raw MeV / 2 / 1000);
                                       // negative flags a relativistic BW
    double charge = 0.0;               // ZPART
    double ctau = 0.0;                 // cm (special flag values: -1, -10, ...)
    int lundId = 0;
    bool decaysInOrder = false;        // numdecays < 0 in the file
    std::vector<DecayMode> decays;
};

struct Config {
    std::string title;
    long nevent = 0;                   // events to generate
    long ngevent = 0;                  // "good" events target
    double pbeamMin = 0.0, pbeamMax = 0.0;  // GeV/c
    int iprule = 0;                    // beam-photon generation rule
    std::string fluxFile;              // external beam flux (iprule=3)
    int nprint = 0;                    // events echoed to terminal
    bool hackFlag = false;
    bool doGluex = false;              // write ASCII event output
    bool doCuts = false;               // call the event analyzer
    std::string prefix;                // output file name prefix
    long maxOutputPerFile = 0;         // events per output file
    double dpOverPOut = 0.0;           // instrumental dp/p smearing
    double pfermi = 0.0;               // GeV/c, quasi-free calculations
    int targetGeom = 0;                // ICODE {0,1,2} = {none,rect,cyl}
    double size1 = 0, size2 = 0, size3 = 0, size4 = 0;  // cm
    int vertexCode = 0;                // JCODE {0,1} = {point,diffuse}
    double xsigma = 0.0, ysigma = 0.0; // cm
    int irestart = 0;                  // 0 = fresh RNG seed
    double reserved = 0.0;

    // Particle table; particles[0] is the beam, particles[1] the target,
    // particles[2] the beam+target pseudo-particle (1-based index i in the
    // Fortran maps to particles[i-1]).
    std::vector<ParticleDef> particles;

    // Find table position by Lund ID (ignoring sign); -1 when absent.
    int findByLund(int lundId) const;
};

// Parse a .def file; throws std::runtime_error with a descriptive message on
// structural errors (mirroring the Fortran's hard exits).
Config parseDefFile(const std::string& path);

}  // namespace mcgen
