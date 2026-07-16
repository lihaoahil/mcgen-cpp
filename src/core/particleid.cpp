#include "core/particleid.h"

namespace mcgen {

int lundToGeant(int lund) {
    switch (lund) {
        case 1: case -1: return 1;    // photon
        case -7: return 2;            // positron
        case 7: return 3;             // electron
        case 2: return 4;             // neutrino
        case -9: return 5;            // mu+
        case 9: return 6;             // mu-
        case 23: case -23: return 7;  // pi0
        case 17: return 8;            // pi+
        case -17: return 9;           // pi-
        case 18: return 11;           // K+
        case -18: return 12;          // K-
        case 42: return 13;           // neutron
        case 41: return 14;           // proton
        case -41: return 15;          // anti-proton
        case 19: case -19: return 16; // K short
        case 24: return 17;           // eta
        case 57: return 18;           // Lambda
        case -57: return 26;          // anti-Lambda
        case 43: return 19;           // Sigma+
        case 44: case -44: return 20; // Sigma0
        case 45: return 21;           // Sigma-
        case 46: return 23;           // Xi-
        case -46: return 31;          // anti-Xi+
        case 47: return 22;           // Xi0
        case -47: return 30;          // anti-Xi0
        case 33: return 57;           // rho0
        case 27: return 58;           // rho+
        case -27: return 59;          // rho-
        case 34: return 60;           // omega
        case 36: return 61;           // eta prime
        case 35: return 62;           // phi
        case 20: case 21: return 0;   // K*0, K*+ (decayed here, not in GEANT)
        case 91: return 991;
        case 92: return 992;
        case 93: return 993;
        case 94: return 994;
        case 95: return 995;
        case 96: return 996;
        case 97: return 997;
        case 98: return 998;
        case 99: return 999;
        case 100: return 9100;
        case 101: case -101: return 9101;
        default: return -1;
    }
}

}  // namespace mcgen
