#include "analysis/rootanalyzer.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

#include "core/dynamics.h"
#include "core/fourvector.h"

namespace mcgen {

namespace {

constexpr double kRade = 57.29577951;
constexpr double kClight = 29.9792458;      // cm/ns
constexpr double kMassProton = 0.938272081;
constexpr double kMassLambda = 1.11568;
constexpr double kMassK0 = 0.497614;

struct Book1 {
    int id;
    const char* title;
    int nx;
    double xlo, xhi;
};
struct Book2 {
    int id;
    const char* title;
    int nx;
    double xlo, xhi;
    int ny;
    double ylo, yhi;
};

// Booking tables transcribed from mc_analyze.F (PAW markup -> TLatex).
const Book1 kBook1[] = {
    {28, "Beam Momentum (GeV/c)", 100, 3.5, 12.0},
    {2800, "Beam Momentum (MeV/c)", 250, 3500., 4500.},
    {131, "P1 Lab Momentum (GeV/c)", 100, 0., 8.0},
    {132, "P2 Lab Momentum (GeV/c)", 100, 0., 8.0},
    {133, "#bar{p} Lab Momentum (GeV/c)", 100, 0., 8.0},
    {135, "p1 Lab Angle (deg)", 100, 0., 50.0},
    {136, "p2 Lab Angle (deg)", 100, 0., 50.0},
    {137, "#bar{p} Lab Angle (deg)", 100, 0., 50.0},
    {138, "Beam (GeV) - 0 FCAL hits", 150, 3.5, 6.5},
    {139, "Beam (GeV) - 1 FCAL hits", 150, 3.5, 6.5},
    {140, "Beam (GeV) - 2 FCAL hits", 150, 3.5, 6.5},
    {141, "Beam (GeV) - 3 FCAL hits", 150, 3.5, 6.5},
    {142, "Beam (GeV)", 150, 3.5, 6.5},
    {148, "FCAL (GeV) - 0 FCAL hits", 300, 0.01, 2.5},
    {149, "FCAL (GeV) - 1 FCAL hits", 300, 0.01, 2.5},
    {150, "FCAL (GeV) - 2 FCAL hits", 300, 0.01, 2.5},
    {151, "FCAL (GeV) - 3 FCAL hits", 300, 0.01, 2.5},
    {152, "FCAL (GeV) - ALL", 300, 0.01, 2.5},
    {153, "CM momentum (GeV/c)", 200, 0.00, 0.8},
    {154, "CM pairwise momentum (GeV/c)", 200, 0.00, 0.75},
    {155, "p1 cos #theta_{LAB}", 50, 0.5, 1.0},
    {156, "p2 cos #theta_{LAB}", 50, 0.5, 1.0},
    {157, "#bar{p} cos #theta_{LAB}", 50, 0.5, 1.0},
    {8611, "Momentum in p#bar{p} pair rest frame (GeV/c)", 100, 0.0, 1.5},
    {8616, "Mom. p#bar{p} (GeV/c) (q-weighted)", 100, 0.0, 1.5},
    {8617, "Mom. p#bar{p} (GeV/c) (q and k weighted)", 100, 0.0, 1.5},
    {8851, "IM (p1 #bar{p}) (GeV/c^{2})", 880, 1.8, 4.0},
    {8852, "IM (p2 #bar{p}) (GeV/c^{2})", 880, 1.8, 4.0},
    {8871, "IM (p1 p2) (GeV/c^{2})", 100, 1.8, 4.0},
    {8857, "p_{p1} (c.m.) (GeV/c)", 100, 0.0, 2.0},
    {8858, "p_{p2} (c.m.) (GeV/c)", 100, 0.0, 2.0},
    {8860, "p_{p1} (c.m.) (RESCALED) (GeV/c)", 100, 0.0, 2.0},
    {8861, "p_{p1} (c.m.) (GeV/c) (RESCALED)", 100, 0.0, 2.0},
    {8862, "p_{p2} (c.m.) (GeV/c) (RESCALED)", 100, 0.0, 2.0},
    {8864, "p_{#bar{p}} (c.m.) (GeV/c)", 100, 0.0, 2.0},
    {8865, "p1 cos #theta_{cm}", 200, -1., 1.0},
    {8866, "p2 cos #theta_{cm}", 200, -1., 1.0},
    {8867, "#bar{p} cos #theta_{cm}", 200, -1., 1.0},
    {865, "p1 cos #theta_{cm}", 200, -1., 1.0},
    {866, "p2 cos #theta_{cm}", 200, -1., 1.0},
    {867, "#bar{p} cos #theta_{cm}", 200, -1., 1.0},
    {1865, "p1 cos #theta_{cm}", 200, -1., 1.0},
    {1866, "p2 cos #theta_{cm}", 200, -1., 1.0},
    {1867, "#bar{p} cos #theta_{cm}", 200, -1., 1.0},
    {2865, "p1 cos #theta_{cm}", 200, -1., 1.0},
    {2866, "p2 cos #theta_{cm}", 200, -1., 1.0},
    {2867, "#bar{p} cos #theta_{cm}", 200, -1., 1.0},
    {3865, "p1 cos #theta_{cm}", 200, -1., 1.0},
    {3866, "p2 cos #theta_{cm}", 200, -1., 1.0},
    {3867, "#bar{p} cos #theta_{cm}", 200, -1., 1.0},
    {8875, "p1 #theta_{cm}", 90, 0., 180.0},
    {8876, "p2 #theta_{cm}", 90, 0., 180.0},
    {8877, "#bar{p} #theta_{cm}", 90, 0., 180.0},
    {8884, "#theta_{cm} IM(p1 pbar)", 90, 0., 180.0},
    {8885, "#theta_{cm} IM(p2 pbar)", 90, 0., 180.0},
    {8886, "cos#theta_{bb} Anti-proton (Adair)", 100, -1.0, 1.0},
    {8887, "cos#theta_{bb} Anti-proton (Hel.)", 100, -1.0, 1.0},
    {8888, "cos#theta_{bb} Anti-proton (GJ)", 100, -1.0, 1.0},
    {8889, "cos #theta_{cm} IM(p1 pbar)", 100, -1.0, 1.0},
    {8890, "cos #theta_{cm} IM(p2 pbar)", 100, -1.0, 1.0},
    {8891, "cos #theta_{cm} p1", 100, -1.0, +1.0},
    {8892, "cos #theta_{cm} p2", 100, -1.0, +1.0},
    {8893, "cos #theta_{cm} pbar", 100, -1.0, +1.0},
    {8895, "-(t) (GeV^{2}) #bar{p} p", 100, 0.0, 5.0},
    {8896, "-(t - t_{min}) (GeV^{2}) #bar{p} p", 100, 0.0, 5.0},
    {8902, "IM (p2 #bar{p}) (GeV/c^{2}) (q-weighted)", 100, 1.8, 4.0},
    {8904, "IM (p2 #bar{p}) (GeV/c^{2}) (q and k weighted)", 100, 1.8, 4.0},
    {8905, "- t_{min} (GeV^{2}) #bar{p} p", 100, 0.0, 5.0},
    {601, "X Vertex Distribution", 100, -3., +3.},
    {602, "Y Vertex Distribution", 100, -3., +3.},
    {603, "Z Vertex Distribution", 700, -50., +20.},
    {604, "X Primary Distribution (cm)", 100, -3., +3.},
    {605, "Y Primary Distribution", 100, -3., +3.},
    {606, "Z Primary Distribution", 100, 0.0, 200.},
    {607, "X #Lambda Distribution (cm)", 100, -3., +3.},
    {608, "Y #Lambda Distribution", 100, -3., +3.},
    {609, "Z #Lambda Distribution", 100, 40., 200.},
    {610, "X #bar{#Lambda} Distribution (cm)", 100, -3., +3.},
    {611, "Y #bar{#Lambda} Distribution", 100, -3., +3.},
    {612, "Z #bar{#Lambda} Distribution", 100, 40., 200.},
    {613, "#Lambda Lab Pathlength (cm)", 100, 0., 100.},
    {614, "#bar{#Lambda} Lab Pathlength (cm)", 100, 0., 100.},
    {615, "#Lambda to #bar{#Lambda} Lab Distance (cm)", 100, 0., 100.},
    {616, "#Lambda Rest Lifetime (ns)", 200, 0., 2.0},
    {617, "#bar{#Lambda} Rest Lifetime (ns)", 200, 0., 2.0},
    {9601, "cos#theta_{Y} of p in #Lambda frm", 50, -1.0, 1.0},
    {9602, "cos#theta_{Y} of #bar{p} in #bar{#Lambda} frm.", 50, -1.0, 1.0},
    {9603, "cos#theta of #Lambda in pair frm.", 50, -1.0, 1.0},
    {9604, "cos#theta of #bar{#Lambda} in pair frm.", 50, -1.0, 1.0},
    {9605, "cos#theta_{X} of p in #Lambda frm.", 50, -1.0, 1.0},
    {9606, "cos#theta_{Y} of p in #Lambda frm.", 50, -1.0, 1.0},
    {9607, "cos#theta_{Z} of p in #Lambda frm.", 50, -1.0, 1.0},
    {9608, "cos#theta_{X} of #bar{p} in #bar{#Lambda} frm.", 50, -1.0, 1.0},
    {9609, "cos#theta_{Y} of #bar{p} in #bar{#Lambda} frm.", 50, -1.0, 1.0},
    {9610, "cos#theta_{Z} of #bar{p} in #bar{#Lambda} frm.", 50, -1.0, 1.0},
    {9611, "Momentum in #Lambda#bar{#Lambda} pair rest frame (GeV/c)", 100,
     0.0, 1.5},
    {9612, "p #bar{#Lambda} in pair rest frame (GeV/c)", 100, 0.0, 1.5},
    {9613, "cos#theta of #Lambda in pair frm.", 50, -1.0, 1.0},
    {9614, "cos#theta of #bar{#Lambda} in pair frm.", 50, -1.0, 1.0},
    {9615, "Momentum in #Lambda#bar{#Lambda} pair rest frame (GeV/c)", 100,
     0.0, 1.5},
    {9616, "Mom. #Lambda#bar{#Lambda} (GeV/c) (q-weighted)", 100, 0.0, 1.5},
    {9617, "Mom. #Lambda#bar{#Lambda} (GeV/c) (q and k weighted)", 100, 0.0,
     1.5},
    {9801, "Lambda IM (p #pi^{-})", 100, 1.05, 1.15},
    {9802, "Anti-Lambda IM (#bar{p} #pi^{+})", 100, 1.05, 1.15},
    {9851, "IM (#bar{#Lambda}#Lambda) (GeV/c^{2})", 100, 2.0, 4.0},
    {9852, "IM (#bar{#Lambda}p) (GeV/c^{2})", 100, 2.0, 4.0},
    {9860, "IM (#Lambda p) (GeV/c^{2})", 100, 2.0, 4.0},
    {9853, "IM (#bar{#Lambda}#Lambda) (GeV/c^{2}) (q-weighted)", 100, 2.0,
     4.0},
    {9858, "IM (#bar{#Lambda}#Lambda) (GeV/c^{2}) (q and k weighted)", 100,
     2.0, 4.0},
    {9856, "-t' (GeV^{2}) #gamma #bar{#Lambda} #Lambda", 100, 0.0, 3.0},
    {9857, "-t' (GeV^{2}) #gamma #bar{#Lambda} p", 100, 0.0, 3.0},
    {9859, "-t' (GeV^{2}) p p recoil", 100, 0.0, 3.0},
    {9861, "-t' (GeV^{2}) #gamma #Lambda", 100, 0.0, 3.0},
    {9862, "-t' (GeV^{2}) #gamma #bar{#Lambda}", 100, 0.0, 3.0},
    {9885, "cos #theta_{cm} Lambda", 100, -1.0, 1.0},
    {9886, "cos #theta_{cm} Anti-Lambda", 100, -1.0, 1.0},
    {9887, "Overall Missing Mass (GeV)", 100, -1.E-1, +5.E-2},
    {9888, "cos #theta_{cm} Between Planes", 100, -1.0, 1.0},
    {9889, "#theta_{cm} Between Planes", 18, 0.0, 180.0},
    {9620, "Overall Missing Energy (GeV)", 100, -0.5, 0.5},
    {9999, "Momentum test", 200, -1.0, 1.0},
    {9901, "Mom p1 (GeV) (lab)", 100, 0.0, 2.0},
    {9902, "Mom p2 (GeV) (lab)", 100, 0.0, 8.0},
    {9903, "Mom p-bar (GeV) (lab)", 100, 0.0, 8.0},
    {9904, "Mom #pi^{+} (GeV) (lab)", 100, 0.0, 2.0},
    {9905, "Mom #pi^{-} (GeV) (lab)", 100, 0.0, 2.0},
    {9906, "cos#theta p1 (lab)", 100, 0.5, 1.0},
    {9907, "cos#theta p2 (lab)", 100, 0.5, 1.0},
    {9908, "cos#theta p-bar (lab)", 100, 0.5, 1.0},
    {9909, "cos#theta #pi^{+} (lab)", 100, 0.5, 1.0},
    {9910, "cos#theta #pi^{-} (lab)", 100, 0.5, 1.0},
    {9911, "Mom p1 (GeV) (cm)", 100, 0.0, 2.0},
    {9912, "Mom p2 (GeV) (cm)", 100, 0.0, 2.0},
    {9913, "Mom p-bar (GeV) (cm)", 100, 0.0, 2.0},
    {9914, "Mom #pi^{+} (GeV) (cm)", 100, 0.0, 1.0},
    {9915, "Mom #pi^{-} (GeV) (cm)", 100, 0.0, 1.0},
    {9916, "cos#theta p1 (cm)", 100, -1., 1.0},
    {9917, "cos#theta p2 (cm)", 100, -1., 1.0},
    {9918, "cos#theta p-bar (cm)", 100, -1., 1.0},
    {9919, "cos#theta #pi^{+} (cm)", 100, -1., 1.0},
    {9920, "cos#theta #pi^{-} (cm)", 100, -1., 1.0},
    {10901, "Mom proton (GeV) (lab)", 100, 0.0, 10.0},
    {10902, "Mom #pi^{+} 1 (GeV) (lab)", 100, 0.0, 6.0},
    {10903, "Mom #pi^{+} 2 (GeV) (lab)", 100, 0.0, 6.0},
    {10904, "Mom #pi^{-} 1 (GeV) (lab)", 100, 0.0, 6.0},
    {10905, "Mom #pi^{-} 2 (GeV) (lab)", 100, 0.0, 6.0},
    {10906, "cos#theta proton (lab)", 100, 0.5, 1.0},
    {10907, "cos#theta #pi^{+} 1 (lab)", 100, 0.5, 1.0},
    {10908, "cos#theta #pi^{+} 2 (lab)", 100, 0.5, 1.0},
    {10909, "cos#theta #pi^{-} 1 (lab)", 100, 0.5, 1.0},
    {10910, "cos#theta #pi^{-} 2 (lab)", 100, 0.5, 1.0},
    {10911, "Mom proton (GeV) (cm)", 100, 0.0, 2.0},
    {10912, "Mom #pi^{+} 1 (GeV) (cm)", 100, 0.0, 2.0},
    {10913, "Mom #pi^{+} 2 (GeV) (cm)", 100, 0.0, 2.0},
    {10914, "Mom #pi^{-} 1 (GeV) (cm)", 100, 0.0, 2.0},
    {10915, "Mom #pi^{-} 2 (GeV) (cm)", 100, 0.0, 2.0},
    {10916, "cos#theta proton (cm)", 100, -1., 1.0},
    {10917, "cos#theta #pi^{+} 1 (cm)", 100, -1., 1.0},
    {10918, "cos#theta #pi^{+} 2 (cm)", 100, -1., 1.0},
    {10919, "cos#theta #pi^{-} 1 (cm)", 100, -1., 1.0},
    {10920, "cos#theta #pi^{-} 2 (cm)", 100, -1., 1.0},
    {10951, "cos#theta K1", 100, -1.0, 1.0},
    {10952, "cos#theta K2", 100, -1.0, 1.0},
    {10885, "cos #theta_{cm} K_{s} 1", 100, -1.0, 1.0},
    {10886, "cos #theta_{cm} K_{s} 2", 100, -1.0, 1.0},
    {10887, "Overall Missing Mass^{2} (GeV)", 100, -2.E-1, +2.E-1},
    {10888, "Overall Missing Mass^{2} (GeV)", 100, -10., +10},
    {10620, "Overall Missing Energy (GeV)", 100, -1.0, 1.0},
    {10630, "IM (K_{s}1 K_{s}2) (GeV/c^{2})", 200, 0.9500, 3.95},
    {10631, "IM (p K_{s}2) (GeV/c^{2})", 200, 1.3, 4.3},
    {10801, "K_{s} 1 IM (#pi^{+} #pi^{-})", 100, 0.300, 0.700},
    {10802, "K_{s} 2 IM (#pi^{+} #pi^{-})", 100, 0.300, 0.700},
    {10851, "IM (K_{s} 1 K_{s} 2) (GeV/c^{2})", 200, 0.95, 3.95},
    {10854, "IM (K_{s} 1 K_{s} 2) (GeV/c^{2})", 1000, 0.95, 3.95},
    {10852, "IM (K_{s} 1 proton) (GeV/c^{2})", 200, 1.3, 4.3},
    {10853, "IM (K_{s} 2 proton) (GeV/c^{2})", 200, 1.3, 4.3},
    {10451, "IM (K_{s} 1 K_{s} 2) (GeV/c^{2})", 200, 0.95, 3.95},
    {10452, "IM (K_{s} 1 proton) (GeV/c^{2})", 200, 1.3, 4.3},
    {10453, "IM (K_{s} 2 proton) (GeV/c^{2})", 200, 1.3, 4.3},
    {1604, "X Primary Distribution (cm)", 100, -3., +3.},
    {1605, "Y Primary Distribution", 100, -3., +3.},
    {1606, "Z Primary Distribution", 100, 0.0, 200.},
    {1607, "X K^{0}_{1} Distribution (cm)", 100, -3., +3.},
    {1608, "Y K^{0}_{1} Distribution", 100, -3., +3.},
    {1609, "Z K^{0}_{1} Distribution", 100, 0.0, 200.},
    {1610, "X K^{0}_{2} Distribution (cm)", 100, -3., +3.},
    {1611, "Y K^{0}_{2} Distribution", 100, -3., +3.},
    {1612, "Z K^{0}_{2} Distribution", 100, 0.0, 200.},
    {1613, "K^{0}_{1} Lab Pathlength (cm)", 100, 0., 100.},
    {1614, "K^{0}_{2} Lab Pathlength (cm)", 100, 0., 100.},
    {1615, "K^{0} Vertex Separation (cm)", 100, 0., 100.},
    {1616, "K^{0}_{1} Rest Lifetime (ns)", 200, 0., 0.5},
    {1617, "K^{0}_{2} Rest Lifetime (ns)", 200, 0., 0.5},
    {1619, "K^{0} Rest Lifetime Sum (ns)", 200, 0., 1.0},
    {1622, "K^{0} Lab Distance Sum (cm)", 100, 0., 100.},
    {1624, "Van Hove Angle (deg)", 180, 0., 360.},
};

const Book2 kBook2[] = {
    {134, "P2 vs. P1 LAB Momenta (GeV/c)", 50, 0., 8.0, 50, 0., 8.0},
    {209, "PROTON #theta (deg) vs. p (GeV/c)", 100, 0., 8., 100, 0., 90.},
    {230, "ANTI-PROTON #theta (deg) vs. p (GeV/c)", 100, 0., 8., 100, 0., 90.},
    {231, "beta vs. p (GeV/c) - ALL", 100, 0., 3.0, 100, 0.4, 1.1},
    {232, "beta vs. p (GeV/c) - 0 FCAL hits", 100, 0., 3.0, 100, 0.4, 1.1},
    {233, "beta vs. p (GeV/c) - 1 FCAL hits", 100, 0., 3.0, 100, 0.4, 1.1},
    {234, "beta vs. p (GeV/c) - 2 FCAL hits", 100, 0., 3.0, 100, 0.4, 1.1},
    {235, "beta vs. p (GeV/c) - 3 FCAL hits", 100, 0., 3.0, 100, 0.4, 1.1},
    {8850, "IM (p2 #bar{p}) vs IM (p1 p2)", 100, 1.6, 4.0, 100, 1.6, 4.0},
    {8853, "IM (p2 #bar{p}) vs IM (p1 #bar{p})", 100, 1.6, 4.0, 100, 1.6, 4.0},
    {8854, "IM^{2} (p2 #bar{p}) vs IM^{2} (p1 #bar{p})", 100, 3.0, 11.0, 100,
     3.0, 11.0},
    {8855, "Longitudinal Momentum (van Hove) Plot", 100, -2.5, 2.5, 100, -2.5,
     2.5},
    {8856, "Proton Momenta P1 vs. P2 in c.m. (GeV/c)", 100, 0.0, 2.0, 100, 0.0,
     2.0},
    {8859, "IM^{2} (p2 #bar{p}) vs IM^{2} (p1 #bar{p})", 50, 3.0, 11.0, 50,
     3.0, 11.0},
    {8863, "IM^{2} (p2 #bar{p}) vs p_{p1} RESCALED", 50, 0.0, 2.0, 50, 3.0,
     11.0},
    {8868, "cos #theta_{cm} (p1 vs. p2)", 50, -1., 1.0, 50, -1., 1.0},
    {8869, "cos #theta_{cm} (pbar vs. p1)", 50, -1., 1.0, 50, -1., 1.0},
    {8870, "cos #theta_{cm} (pbar vs. p2)", 50, -1., 1.0, 50, -1., 1.0},
    {8878, "#theta_{cm} (p1 vs. p2)", 90, 0., 180.0, 90, 0., 180.0},
    {8879, "#theta_{cm} (pbar vs. p1)", 90, 0., 180.0, 90, 0., 180.0},
    {8880, "#theta_{cm} (pbar vs. p2)", 90, 0., 180.0, 90, 0., 180.0},
    {8881, "#theta_{cm} (deg) vs. p_{cm} Proton1", 100, 0.0, 2.0, 90, 0, 180.},
    {8882, "#theta_{cm} (deg) vs. p_{cm} Proton2", 100, 0.0, 2.0, 90, 0, 180.},
    {8883, "#theta_{cm} (deg) vs. p_{cm} Antiproton", 100, 0.0, 2.0, 90, 0,
     180.},
    {8894, "cos#theta_{cm} (p2 vs. p1) BEFORE SWAPPING", 50, -1.0, 1.0, 50,
     -1.0, 1.0},
    {8897, "Beam Energy (GeV) vs. IM(p_{2} #bar{p})", 220, 1.8, 4.0, 90, 3.0,
     12.0},
    {8898, "-(t - t_{min}) (GeV^{2}) #bar{p} p vs. IM(p_{2} #bar{p})", 110,
     1.8, 4.0, 50, 0.0, 3.0},
    {8899, "-t (GeV^{2}) #bar{p} p vs. IM(p_{2} #bar{p})", 110, 1.8, 4.0, 50,
     0.0, 3.0},
    {8900, "-t (GeV^{2}) vs. #bar{p} cos #theta_{cm}", 50, -1., +1.0, 50, 0.0,
     3.0},
    {8901, "Beam Energy (GeV) vs. #bar{p} cos #theta_{cm}", 50, -1., +1.0, 90,
     3.0, 12.0},
    {618, "#Lambda vs. #bar{#Lambda} Rest Lifetimes (ns)", 30, 0., 1.0, 30, 0.,
     1.0},
    {9850, "IM (p #bar{#Lambda}) vs IM (p #Lambda)", 100, 2.0, 4.0, 100, 2.0,
     4.0},
    {9855, "IM (p #bar{#Lambda}) vs IM (#bar{#Lambda} #Lambda)", 100, 2.0, 4.0,
     100, 2.0, 4.0},
    {9849, "IM^{2} (p #bar{#Lambda}) vs IM^{2} (#bar{#Lambda} #Lambda) "
           "(Dalitz Plot)",
     100, 4.5, 10.5, 100, 4.0, 10.0},
    {9854, "IM^{2} (p #bar{#Lambda}) vs IM^{2} (p #Lambda) (Dalitz Plot)", 100,
     4.0, 10.0, 100, 4.0, 10.0},
    {9869, "cos #theta_{cm} anti-Lambda vs. Lambda", 50, -1., 1.0, 50, -1.,
     1.0},
    {9870, "Long. Mom. (van Hove)", 100, -2.5, 2.5, 100, -2.5, 2.5},
    {9871, "Long. Mom. Inv. Mass (Schu Plot)", 100, -3.5, 3.5, 100, -3.5, 3.5},
    {9872, "Long. Mom. Kin. Enrgy. (Schu2 Plot)", 100, -2.5, 2.5, 100, -2.5,
     2.5},
    {9873, "Long. Mom. Trans. Mom. (Schu3 Plot)", 100, -2.5, 2.5, 100, -2.5,
     2.5},
    {9874, "I.M. (by sector pair) vs. #omega (degrees)", 100, 0.0, 360, 50,
     2.0, 3.5},
    {9880, "p_{cm} (GeV/c) anti-Lambda vs. Lambda", 100, 0.0, 2.0, 100, 0, 2.},
    {9881, "cos #theta_{cm} vs. p_{cm} Lambda", 100, 0.0, 2.0, 100, -1.0,
     +1.0},
    {9883, "cos #theta_{cm} vs. p_{cm} anti-Lambda", 100, 0.0, 2.0, 100, -1.0,
     +1.0},
    {9884, "cos #theta_{cm} vs. p_{cm} (#bar{#Lambda} #Lambda Pair)", 100, 0.0,
     2.0, 100, -1.0, +1.0},
    {9621, "Missing Energy vs. Missing Mass", 50, -1.E-1, +5.E-2, 50, -0.5,
     0.5},
    {9625, "-t' (GeV^{2}) #Lambda vs. IM #Lambda#bar{#Lambda}", 100, 2.2, 4.0,
     100, 0.0, 3.0},
    {9626, "-t' (GeV^{2}) #bar{#Lambda} vs. IM #Lambda#bar{#Lambda}", 100, 2.2,
     4.0, 100, 0.0, 3.0},
    {9627, "-t' (GeV^{2}) p p vs. IM #Lambda#bar{#Lambda}", 100, 2.2, 4.0, 100,
     0.0, 3.0},
    {9628, "-t' (GeV^{2}) #Lambda#bar{#Lambda} vs. IM #Lambda#bar{#Lambda}",
     100, 2.2, 4.0, 100, 0.0, 3.0},
    {9630, "-t (GeV^{2}) p p vs. -t (GeV^{2}) #gamma #Lambda #Lambda", 100,
     -10.0, 10.0, 100, -10.0, 10.0},
    {9631, "-t' (GeV^{2}) p p vs. -t' (GeV^{2}) #gamma #Lambda #Lambda", 100,
     -10.0, 10.0, 100, -10.0, 10.0},
    {9890, "Correlation Matrix #bar{#Lambda} vs #Lambda", 5, 0.5, 5.5, 3, 0.5,
     3.5},
    {10921, "t_{K1} (GeV^{2}) vs. IM(K_{2} p) (GeV)", 100, 0.95, 3.95, 100,
     0.0, 1.5},
    {10922, "t_{K1}, t_{K2} (GeV^{2}) vs. IM(K_{1} K_{2}) (GeV)", 100, 0.95,
     3.95, 100, 0.0, 7.5},
    {10923, "t_{K1 K2} (GeV^{2}) vs. IM(K_{2} p) (GeV)", 100, 0.95, 3.95, 100,
     0.0, 4.0},
    {10924, "t_{K1 K2} (GeV^{2}) vs. IM(K_{1} K_{2}) (GeV)", 100, 0.95, 3.95,
     100, 0.0, 1.5},
    {10931, "cos #theta_{K1} vs. IM(K_{2} p) (GeV)", 100, 0.95, 3.95, 100,
     -1.0, +1.0},
    {10932, "cos#theta_{K1} cos#theta_{K2} vs. IM(K_{1} K_{2}) (GeV)", 100,
     0.95, 3.95, 100, -1.0, +1.0},
    {10933, "cos#theta_{K1 K2} vs. IM(K_{2} p) (GeV)", 100, 0.95, 3.95, 100,
     -1.0, +1.0},
    {10934, "cos#theta_{K1 K2} vs. IM(K_{1} K_{2}) (GeV)", 100, 0.95, 3.95,
     100, -1.0, +1.0},
    {10941, "-t K2 (GeV^{2}) vs. -t K1 (GeV^{2})", 100, 0.0, 5.0, 100, 0.0,
     5.0},
    {10942, "cos#theta K2 vs. cos#theta K1", 100, -1.0, 1.0, 100, -1.0, +1.0},
    {10889, "cos#theta K^{0}_{2} vs K^{0}_{1} (cm)", 100, -1.0, 1.0, 100, -1.0,
     1.0},
    {10621, "Missing Energy vs. Missing Mass", 100, -5.E-1, +5.E-1, 100, -1.5,
     1.5},
    {10622, "Mass K_{s} 2 vs. Mass K_{s} 1", 100, 0.3, 0.7, 100, 0.3, 0.7},
    {10623, "IM (p K_{s}2) vs. IM(p K_{s}1)", 100, 1.3, 4.3, 100, 1.3, 4.3},
    {10624, "IM (p K_{s}2) vs. IM(K_{s}1 K_{s}2)", 100, 0.4, 4.4, 100, 1.3,
     4.3},
    {1618, "K^{0}_{1} vs. K^{0}_{2} Rest Lifetimes (ns)", 30, 0., 0.5, 30, 0.,
     0.5},
    {1620, "Longitudinal Momentum (van Hove) Plot", 100, -2.5, 2.5, 100, -2.5,
     2.5},
    {1621, "van Hove Angle (deg) vs. IM (Ks p)", 100, 1.3, 4.3, 100, 0.0,
     360.},
    {1623, "Rest Lifetimes vs. Lab Distance Sums", 30, 0., 10.0, 30, 0., 0.1},
};

// GETANGLES: scalar momentum, polar angle (deg), cos(theta), phi (deg).
void getangles(const FourVec& p, double& ppart, double& theta,
               double& costheta, double& phi) {
    ppart = p.p();
    if (ppart > 0.0) {
        costheta = std::clamp(p.pz / ppart, -1.0, 1.0);
        theta = std::acos(costheta) * kRade;
        phi = std::atan2(p.py, p.px) * kRade;
        if (phi < 0.0) phi += 360.0;
    } else {
        costheta = 2.0;
        theta = 1000.0;
        phi = 1000.0;
    }
}

// MISSINGMASS: missing mass and its square of p1 - p2.
void missingmass(const FourVec& p1, const FourVec& p2, double& rmiss,
                 double& rmiss2) {
    const double px = p1.px - p2.px, py = p1.py - p2.py, pz = p1.pz - p2.pz;
    const double ee = p1.e - p2.e;
    rmiss2 = ee * ee - (px * px + py * py + pz * pz);
    rmiss = rmiss2 > 0.0 ? std::sqrt(rmiss2) : -0.000001;
}

// GETDECAYANGLE
void getdecayangle(const FourVec& p, const Vec3& beta, double& costh,
                   double& theta) {
    const double bmag = std::sqrt(beta[0] * beta[0] + beta[1] * beta[1] +
                                  beta[2] * beta[2]);
    costh = (p.px * beta[0] + p.py * beta[1] + p.pz * beta[2]) /
            (p.p() * bmag);
    theta = std::acos(costh) * kRade;
}

Vec3 getbeta(const FourVec& p) {
    return {p.px / p.e, p.py / p.e, p.pz / p.e};
}

// Van Hove sector from the three longitudinal momenta (0 = unsortable).
int vanHoveSector(double q1, double q2, double q3) {
    if (q1 >= 0 && q2 <= 0 && q3 >= 0) return 1;
    if (q1 >= 0 && q2 <= 0 && q3 <= 0) return 2;
    if (q1 >= 0 && q2 >= 0 && q3 <= 0) return 3;
    if (q1 <= 0 && q2 >= 0 && q3 <= 0) return 4;
    if (q1 <= 0 && q2 >= 0 && q3 >= 0) return 5;
    if (q1 <= 0 && q2 <= 0 && q3 >= 0) return 6;
    return 0;
}

// IACCEPT_GLUEX: rough GlueX single-track acceptance.  Returns the accept
// flag; acceptweight gets the decay-probability weight.
bool iacceptGluex(double theta, double pmomMeV, int icharge, double bgctau,
                  int lund, double& acceptweight) {
    bool iacc = true;
    acceptweight = std::exp(-100.0 / bgctau);  // drift-chamber flight path
    const double pmomgev = pmomMeV / 1000.0;
    if (icharge != 0) {
        double pmin = 0.15;                       // generic charged
        if (std::abs(lund) == 17) pmin = 0.15;    // pions
        else if (std::abs(lund) == 41) pmin = 0.35;  // (anti)protons
        if (pmomgev < pmin) {
            iacc = false;
            acceptweight = 0.0;
        }
        if (theta > 150.0 || theta < 2.0) {
            iacc = false;
            acceptweight = 0.0;
        }
    } else if (lund == 1) {  // photon
        if (theta < 2.0 || theta > 125.0 || pmomgev < 0.03) {
            iacc = false;
            acceptweight = 0.0;
        }
    }
    return iacc;
}

}  // namespace

struct RootAnalyzer::Impl {
    explicit Impl(const std::string& histFile) {
        std::printf(" First call to event analyzer (C++/ROOT port)\n");
        file = std::make_unique<TFile>(histFile.c_str(), "RECREATE");
        static const char* dirs[3] = {"RAWEVENT", "ACCEPTED", "WEIGHTED"};
        for (int d = 0; d < 3; ++d) {
            TDirectory* dir = file->mkdir(dirs[d]);
            dir->cd();
            for (const auto& b : kBook1) {
                auto* h = new TH1D(Form("h%d", b.id), b.title, b.nx, b.xlo,
                                   b.xhi);
                h1[d][b.id] = h;
            }
            for (const auto& b : kBook2) {
                auto* h = new TH2D(Form("h%d", b.id), b.title, b.nx, b.xlo,
                                   b.xhi, b.ny, b.ylo, b.yhi);
                h2[d][b.id] = h;
            }
        }
        for (auto& row : corrsum) row.fill(0.0);
        for (auto& row : correrr) row.fill(0.0);
    }

    void hf1(int id, double x, double w) {
        auto it = h1[dir].find(id);
        if (it != h1[dir].end()) it->second->Fill(x, w);
    }
    void hf2(int id, double x, double y, double w) {
        auto it = h2[dir].find(id);
        if (it != h2[dir].end()) it->second->Fill(x, y, w);
    }

    std::unique_ptr<TFile> file;
    int dir = 0;  // 0 RAWEVENT, 1 ACCEPTED, 2 WEIGHTED
    std::map<int, TH1D*> h1[3];
    std::map<int, TH2D*> h2[3];

    long icount = 0, ikeep = 0;
    bool firstProcess = true;
    int imechanism = 0;  // Lambda anti-Lambda mechanism, set on iloop 1
    Rng rng{424242};     // for the KsKs randomized swap (RAN in the Fortran)
    // corrsum/correrr(3,5): columns 0-2 = correlation matrix,
    // 3 = Lambda polarization sums, 4 = anti-Lambda polarization sums.
    std::array<std::array<double, 5>, 3> corrsum{}, correrr{};
};

RootAnalyzer::RootAnalyzer(const std::string& histFile)
    : impl_(std::make_unique<Impl>(histFile)) {}

RootAnalyzer::~RootAnalyzer() = default;

bool RootAnalyzer::analyze(const Event& ev, double ppbeam) {
    auto& im = *impl_;
    ++im.icount;
    const int imark = static_cast<int>(ev.size());  // last Fortran index

    // Per-track acceptance.
    std::vector<double> acceptlist(ev.size(), 0.0);
    int nProt = 0, nAntip = 0, nPip = 0, nPim = 0;
    for (size_t i = 1; i < ev.size(); ++i) {
        const auto& q = ev[i];
        const double pmomMeV = q.pDecay.p() * 1000.0;
        double theta, costh, phi, pmag;
        getangles(q.pDecay, pmag, theta, costh, phi);
        const double betamag = q.pDecay.p() / q.pDecay.e;
        double gamma = 1.0 - betamag * betamag;
        gamma = gamma > 0.0 ? 1.0 / std::sqrt(gamma) : 1.e32;
        const double bgctau = betamag * gamma * q.ctau;

        double w = 0.0;
        const bool iacc =
            iacceptGluex(theta, pmomMeV, static_cast<int>(q.charge), bgctau,
                         q.lund, w);
        acceptlist[i] = iacc ? w : 0.0;
        if (iacc) {
            if (q.lund == 41) ++nProt;
            else if (q.lund == -41) ++nAntip;
            else if (q.lund == 17) ++nPip;
            else if (q.lund == -17) ++nPim;
        }
    }

    // Software trigger by event topology (hard-coded, as in the Fortran).
    int itrigmask = 1;
    if (imark == 4 || imark == 5) itrigmask = 2;        // p pbar p
    else if (imark == 8 || imark == 9) itrigmask = 3;   // p pbar p pi+ pi-
    else if (imark == 6 || imark == 11) itrigmask = 4;  // Ks Ks p

    bool iprocess = false;
    if (itrigmask == 1) iprocess = true;
    else if (itrigmask == 2 && nProt >= 2 && nAntip >= 1) iprocess = true;
    else if (itrigmask == 3 && nProt >= 2 && nAntip >= 1 && nPip >= 1 &&
             nPim >= 1)
        iprocess = true;
    else if (nProt >= 1 && nPip >= 2 && nPim >= 2)
        iprocess = true;  // Ks Ks trigger, set after the fact

    const bool iwinner = iprocess;
    if (iwinner && im.firstProcess) {
        std::printf(" The first event to pass the trigger selection\n");
        im.firstProcess = false;
    }

    FourVec pbeamlab{0.0, 0.0, ppbeam, 0.0, ppbeam};

    // Three passes: RAWEVENT (weight 1), ACCEPTED (weight 1, only winners),
    // WEIGHTED (acceptance-weighted; only the KsKs block computes one).
    for (int iloop = 1; iloop <= 3; ++iloop) {
        double weight = 1.0;
        if (iloop == 1) {
            im.dir = 0;
        } else if (iloop == 2) {
            if (!iwinner) return false;
            im.dir = 1;
        } else {
            im.dir = 2;
            weight = 0.0;  // itrigmask 10/11 paths are never set (legacy)
        }

        const double ww = ev[0].p.m;
        const Vec3 beta = getbeta(ev[0].p);

        // ---------------- p p pbar block (3 tracks) ----------------
        if (imark <= 5) {
            for (size_t i = 1; i < ev.size(); ++i) {
                if (ev[i].lund != 41) continue;
                for (size_t j = i + 1; j < ev.size(); ++j) {
                    if (ev[j].lund != 41) continue;
                    for (size_t k = 1; k < ev.size(); ++k) {
                        if (ev[k].lund != -41) continue;
                        if (iloop == 1 && iprocess) ++im.ikeep;

                        // beta = +p/E of the total system boosts into the CM.
                        FourVec p1cm = boost(ev[i].p, beta);
                        FourVec p2cm = boost(ev[j].p, beta);
                        FourVec apcm = boost(ev[k].p, beta);

                        double pm1cm, th1cm, cth1cm, ph1cm;
                        double pm2cm, th2cm, cth2cm, ph2cm;
                        double pmacm, thacm, cthacm, phacm;
                        getangles(p1cm, pm1cm, th1cm, cth1cm, ph1cm);
                        getangles(p2cm, pm2cm, th2cm, cth2cm, ph2cm);
                        getangles(apcm, pmacm, thacm, cthacm, phacm);
                        im.hf2(8894, cth1cm, cth2cm, weight);

                        FourVec pimp1p2, pimp1ap, pimp2ap;
                        const double rmp1p2 = invarMass(ev[i].p, ev[j].p, pimp1p2);
                        const double rmp1ap = invarMass(ev[i].p, ev[k].p, pimp1ap);
                        const double rmp2ap = invarMass(ev[j].p, ev[k].p, pimp2ap);
                        im.hf1(8851, rmp1ap, weight);
                        im.hf1(8852, rmp2ap, weight);
                        im.hf1(8871, rmp1p2, weight);
                        im.hf2(8850, rmp1p2, rmp2ap, weight);
                        im.hf2(8853, rmp1ap, rmp2ap, weight);
                        im.hf2(8854, rmp1ap * rmp1ap, rmp2ap * rmp2ap, weight);
                        im.hf2(8897, rmp2ap, ppbeam, weight);

                        im.hf1(28, ppbeam, weight);
                        im.hf1(2800, ppbeam * 1000, weight);
                        double pv1, th1, cth1, ph1, pv2, th2, cth2, ph2, pva,
                            tha, ctha, pha;
                        getangles(ev[i].p, pv1, th1, cth1, ph1);
                        getangles(ev[j].p, pv2, th2, cth2, ph2);
                        getangles(ev[k].p, pva, tha, ctha, pha);
                        const double beta1 = pv1 / ev[i].p.e;
                        const double beta2 = pv2 / ev[j].p.e;
                        const double betaa = pva / ev[k].p.e;
                        const double ep1 = ev[i].p.e - kMassProton;
                        const double ep2 = ev[j].p.e - kMassProton;
                        const double epa = ev[k].p.e - kMassProton;
                        im.hf1(131, pv1, weight);
                        im.hf1(132, pv2, weight);
                        im.hf1(133, pva, weight);
                        im.hf2(134, pv1, pv2, weight);
                        im.hf1(135, th1, weight);
                        im.hf1(136, th2, weight);
                        im.hf1(137, tha, weight);
                        im.hf1(155, cth1, weight);
                        im.hf1(156, cth2, weight);
                        im.hf1(157, ctha, weight);
                        im.hf2(209, pv1, th1, weight);
                        im.hf2(230, pva, tha, weight);

                        // FCAL hit-count estimate.
                        const double plimit = 4.0, angval = 12.0;
                        int numfcal = 0;
                        if (th1 <= angval) ++numfcal;
                        if (th2 <= angval) ++numfcal;
                        if (tha <= angval) ++numfcal;
                        const int idbeam = 138 + numfcal;
                        const int idcth = numfcal * 1000 + 865;
                        const int idbeta = 232 + numfcal;
                        im.hf1(idbeam, ppbeam, weight);
                        im.hf1(idcth, cth1cm, weight);
                        im.hf1(idcth + 1, cth2cm, weight);
                        im.hf1(idcth + 2, cthacm, weight);
                        double esum = 0.0;
                        if (ppbeam <= plimit) {
                            if (th1 < 10.0) { im.hf2(idbeta, pv1, beta1, weight); im.hf2(231, pv1, beta1, weight); esum += ep1; }
                            if (th2 < 10.0) { im.hf2(idbeta, pv2, beta2, weight); im.hf2(231, pv2, beta2, weight); esum += ep2; }
                            if (tha < 10.0) { im.hf2(idbeta, pva, betaa, weight); im.hf2(231, pva, betaa, weight); esum += epa; }
                        }
                        im.hf1(142, ppbeam, weight);
                        // esumN is nonzero only for the matching hit count,
                        // but all five spectra get an entry (as in HBOOK).
                        for (int nf = 0; nf <= 3; ++nf)
                            im.hf1(148 + nf, nf == numfcal ? esum : 0.0,
                                   weight);
                        im.hf1(152, esum, weight);

                        // CM momenta and pairwise momenta.
                        im.hf1(153, pm1cm, weight);
                        im.hf1(153, pm2cm, weight);
                        im.hf1(153, pmacm, weight);
                        for (auto [a, b] : {std::pair{&p1cm, &p2cm},
                                            {&p1cm, &apcm},
                                            {&p2cm, &apcm}}) {
                            FourVec pair;
                            pair.px = a->px - b->px;
                            pair.py = a->py - b->py;
                            pair.pz = a->pz - b->pz;
                            double pm, th, cth, ph;
                            getangles(pair, pm, th, cth, ph);
                            im.hf1(154, pm, weight);
                        }

                        // t of the beam to the (p2, pbar) pair, in the CM.
                        FourVec pbeamcm = boost(pbeamlab, beta);
                        FourVec pimp1apcm = boost(pimp1ap, beta);
                        FourVec pimp2apcm = boost(pimp2ap, beta);
                        double pmomim1, thim1, cthim1, phim1;
                        double pmomim2, thim2, cthim2, phim2;
                        getangles(pimp1apcm, pmomim1, thim1, cthim1, phim1);
                        getangles(pimp2apcm, pmomim2, thim2, cthim2, phim2);
                        double tt, ttmin, ttprime;
                        gett(pbeamcm, pimp2apcm, tt, ttmin, ttprime);
                        im.hf1(8895, -tt, weight);
                        im.hf1(8896, -ttprime, weight);
                        im.hf1(8905, -ttmin, weight);
                        im.hf2(8898, rmp2ap, -ttprime, weight);
                        im.hf2(8899, rmp2ap, -tt, weight);
                        im.hf2(8900, cthacm, -tt, weight);
                        im.hf2(8901, cthacm, ppbeam, weight);

                        im.hf2(8856, pm2cm, pm1cm, weight);
                        im.hf1(8857, pm1cm, weight);
                        im.hf1(8858, pm2cm, weight);
                        im.hf1(8864, pmacm, weight);
                        im.hf1(8865, cth1cm, weight);
                        im.hf1(8866, cth2cm, weight);
                        im.hf1(8867, cthacm, weight);
                        im.hf2(8868, cth2cm, cth1cm, weight);
                        im.hf2(8869, cth1cm, cthacm, weight);
                        im.hf2(8870, cth2cm, cthacm, weight);
                        im.hf1(8875, th1cm, weight);
                        im.hf1(8876, th2cm, weight);
                        im.hf1(8877, thacm, weight);
                        im.hf2(8878, th2cm, th1cm, weight);
                        im.hf2(8879, th1cm, thacm, weight);
                        im.hf2(8880, th2cm, thacm, weight);
                        im.hf2(8881, pm1cm, th1cm, weight);
                        im.hf2(8882, pm2cm, th2cm, weight);
                        im.hf2(8883, pmacm, thacm, weight);
                        im.hf1(8884, thim1, weight);
                        im.hf1(8885, thim2, weight);
                        im.hf1(8889, cthim1, weight);
                        im.hf1(8890, cthim2, weight);
                        im.hf1(8891, cth1cm, weight);
                        im.hf1(8892, cth2cm, weight);
                        im.hf1(8893, cthacm, weight);

                        // Rescaled Dalitz plot (Egamma = 8.5 -> W = 4.08).
                        const double wtarget = 4.08;
                        const double twomp2 = (2.0 * kMassProton) * (2.0 * kMassProton);
                        const double scale =
                            ((wtarget - kMassProton) * (wtarget - kMassProton) - twomp2) /
                            ((ww - kMassProton) * (ww - kMassProton) - twomp2);
                        const double rmp1ap2hack = (rmp1ap * rmp1ap - twomp2) * scale + twomp2;
                        const double rmp2ap2hack = (rmp2ap * rmp2ap - twomp2) * scale + twomp2;
                        im.hf2(8859, rmp1ap2hack, rmp2ap2hack, weight);
                        const double pcm = (1.0 / (2.0 * ww)) *
                                           std::sqrt(ww * ww * (ww * ww - twomp2));
                        const double ptargetv = (1.0 / (2.0 * wtarget)) *
                                                std::sqrt(wtarget * wtarget *
                                                          (wtarget * wtarget - twomp2));
                        const double sratio = ptargetv / pcm;
                        im.hf1(8861, pm1cm * sratio, weight);
                        im.hf1(8862, pm2cm * sratio, weight);
                        im.hf2(8863, pm1cm * sratio, rmp2ap2hack, weight);

                        // Baryon-baryon rest frame angles.
                        const Vec3 betabb = getbeta(pimp2apcm);
                        FourVec p2bb = boost(p2cm, betabb);
                        FourVec apbb = boost(apcm, betabb);
                        FourVec pbeambb = boost(pbeamcm, betabb);
                        double p2mombb, thp2bb, cthp2bb, php2bb;
                        double apmombb, thapbb, cthapbb, phapbb;
                        getangles(p2bb, p2mombb, thp2bb, cthp2bb, php2bb);
                        getangles(apbb, apmombb, thapbb, cthapbb, phapbb);
                        const double bmag =
                            std::sqrt(betabb[0] * betabb[0] + betabb[1] * betabb[1] +
                                      betabb[2] * betabb[2]);
                        const double cthhel =
                            (apbb.px * betabb[0] + apbb.py * betabb[1] +
                             apbb.pz * betabb[2]) / (apbb.p() * bmag);
                        const double cthgj = dot3(apbb, pbeambb) /
                                             (apbb.p() * pbeambb.p());
                        im.hf1(8886, cthapbb, weight);
                        im.hf1(8887, cthhel, weight);
                        im.hf1(8888, cthgj, weight);

                        // Phase-space-deweighted pair spectra.
                        const double pscale = 0.5;
                        const double reweight = weight / p2mombb * pscale;
                        const double reweight2 = reweight / pmomim2;
                        im.hf1(8611, p2mombb, weight);
                        im.hf1(8616, p2mombb, reweight);
                        im.hf1(8617, p2mombb, reweight2);
                        im.hf1(8902, rmp2ap, reweight);
                        im.hf1(8904, rmp2ap, reweight2);

                        // Van Hove plots, both proton assignments.
                        for (int icombo = 0; icombo < 2; ++icombo) {
                            const double q1 = icombo == 0 ? p1cm.pz : p2cm.pz;
                            const double q2 = apcm.pz;
                            const double q3 = icombo == 0 ? p2cm.pz : p1cm.pz;
                            const int isector = vanHoveSector(q1, q2, q3);
                            if (isector > 0) {
                                const double qq =
                                    0.81650 * std::sqrt(q1 * q1 + q2 * q2 + q3 * q3);
                                const double omega =
                                    std::atan2(-std::sqrt(3.0) * q1,
                                               q1 + 2.0 * q2) + M_PI;
                                im.hf2(9870, qq * std::cos(omega),
                                       qq * std::sin(omega), weight);
                            } else if (iloop == 1) {
                                std::printf(" Cannot sort right... %g %g %g\n",
                                            q1, q2, q3);
                            }
                        }
                    }
                }
            }
        }

        // ------------- Lambda anti-Lambda block (5 tracks) -------------
        // A false return = the Fortran "goto 1000": skip the rest of this
        // iloop pass and move on to the next one.
        if (!lambdaBlock(ev, ppbeam, weight, iloop, iprocess)) continue;
        if (!ksksBlock(ev, ppbeam, weight, iloop, iprocess, acceptlist))
            continue;
    }
    return iwinner;
}

bool RootAnalyzer::lambdaBlock(const Event& ev, double ppbeam, double weight,
                               int iloop, bool iprocess) {
    auto& im = *impl_;
    const int imark = static_cast<int>(ev.size());
    const double ww = ev[0].p.m;
    const Vec3 beta = getbeta(ev[0].p);
    FourVec pbeamlab{0.0, 0.0, ppbeam, 0.0, ppbeam};
    FourVec ptargetlab{0.0, 0.0, 0.0, kMassProton, kMassProton};

    // Mechanism from the fixed particle ordering in the .def (iloop 1 only).
    if (iloop == 1) {
        const int ipart2 = imark >= 2 ? ev[1].lund : 0;
        const int ipart3 = imark >= 3 ? ev[2].lund : 0;
        if (ipart2 == 41) im.imechanism = 1;
        else if (ipart2 == 57 && ipart3 == 96) im.imechanism = 2;
        else if (ipart2 == 96) im.imechanism = 3;
        else if (ipart2 == 94) im.imechanism = 4;
        else if (ipart2 == 57 && ipart3 == 41) im.imechanism = 5;
    }
    const int mech = im.imechanism;

    for (size_t i = 1; i < ev.size(); ++i) {
        if (ev[i].lund != 41) continue;
        for (size_t j = i + 1; j < ev.size(); ++j) {
            if (ev[j].lund != 41) continue;
            for (size_t k = 1; k < ev.size(); ++k) {
                if (ev[k].lund != -41) continue;
                for (size_t l = 1; l < ev.size(); ++l) {
                    if (ev[l].lund != 17) continue;
                    for (size_t m = 1; m < ev.size(); ++m) {
                        if (ev[m].lund != -17) continue;

    if (iloop == 1 && iprocess) ++im.ikeep;

    // Lab momenta / angles.
    double p1momlab, thp1lab, cthp1lab, php1lab;
    double p2momlab, thp2lab, cthp2lab, php2lab;
    double apmomlab, thaplab, cthaplab, phaplab;
    double pipmomlab, thpiplab, cthpiplab, phpiplab;
    double pimmomlab, thpimlab, cthpimlab, phpimlab;
    getangles(ev[i].p, p1momlab, thp1lab, cthp1lab, php1lab);
    getangles(ev[j].p, p2momlab, thp2lab, cthp2lab, php2lab);
    getangles(ev[k].p, apmomlab, thaplab, cthaplab, phaplab);
    getangles(ev[l].p, pipmomlab, thpiplab, cthpiplab, phpiplab);
    getangles(ev[m].p, pimmomlab, thpimlab, cthpimlab, phpimlab);
    im.hf1(9901, p1momlab, weight);
    im.hf1(9902, p2momlab, weight);
    im.hf1(9903, apmomlab, weight);
    im.hf1(9904, pipmomlab, weight);
    im.hf1(9905, pimmomlab, weight);
    im.hf1(9906, cthp1lab, weight);
    im.hf1(9907, cthp2lab, weight);
    im.hf1(9908, cthaplab, weight);
    im.hf1(9909, cthpiplab, weight);
    im.hf1(9910, cthpimlab, weight);

    // Boost to the CM: generated (iloop 1) vs smeared (iloop 2/3) momenta.
    const bool gen = iloop == 1;
    auto pick = [&](size_t idx) { return gen ? ev[idx].p : ev[idx].pDecay; };
    FourVec p1cm = boost(pick(i), beta);
    FourVec p2cm = boost(pick(j), beta);
    FourVec apcm = boost(pick(k), beta);
    FourVec pipcm = boost(pick(l), beta);
    FourVec pimcm = boost(pick(m), beta);
    FourVec pbeamcm = boost(pbeamlab, beta);
    FourVec aplab = pick(k);
    FourVec piplab = pick(l);
    FourVec psystem = boost(ev[0].p, beta);
    FourVec ptargetcm = boost(ptargetlab, beta);

    double p1momcm, thp1cm, cthp1cm, php1cm;
    double p2momcm, thp2cm, cthp2cm, php2cm;
    double apmomcm, thapcm, cthapcm, phapcm;
    double pipmomcm, thpipcm, cthpipcm, phpipcm;
    double pimmomcm, thpimcm, cthpimcm, phpimcm;
    getangles(p1cm, p1momcm, thp1cm, cthp1cm, php1cm);
    getangles(p2cm, p2momcm, thp2cm, cthp2cm, php2cm);
    getangles(apcm, apmomcm, thapcm, cthapcm, phapcm);
    getangles(pipcm, pipmomcm, thpipcm, cthpipcm, phpipcm);
    getangles(pimcm, pimmomcm, thpimcm, cthpimcm, phpimcm);
    im.hf1(9911, p1momcm, weight);
    im.hf1(9912, p2momcm, weight);
    im.hf1(9913, apmomcm, weight);
    im.hf1(9914, pipmomcm, weight);
    im.hf1(9915, pimmomcm, weight);
    im.hf1(9916, cthp1cm, weight);
    im.hf1(9917, cthp2cm, weight);
    im.hf1(9918, cthapcm, weight);
    im.hf1(9919, cthpipcm, weight);
    im.hf1(9920, cthpimcm, weight);

    // Invariant-mass combinations; the pairing depends on the mechanism
    // (the .def must create particles in a fixed order for this to work).
    FourVec alamcm, alamlab;
    const double rmalam = invarMass(apcm, pipcm, alamcm);
    invarMass(aplab, piplab, alamlab);

    Vec3 rvertex, dvertex;
    FourVec plamlab, plamcm, ppplam, ppalam;
    double rmplam, rmpplam, rmpalam;
    const size_t recoilIdx = (mech == 2) ? j : i;
    const size_t decayIdx = (mech == 2) ? i : j;
    rvertex = ev[recoilIdx].rCreate;
    dvertex = ev[decayIdx].rCreate;
    invarMass(pick(decayIdx), pick(m), plamlab);
    if (mech == 2) {
        rmplam = invarMass(p1cm, pimcm, plamcm);
        rmpplam = invarMass(p2cm, plamcm, ppplam);
        rmpalam = invarMass(p2cm, alamcm, ppalam);
    } else {
        rmplam = invarMass(p2cm, pimcm, plamcm);
        rmpplam = invarMass(p1cm, plamcm, ppplam);
        rmpalam = invarMass(p1cm, alamcm, ppalam);
    }

    FourVec plamalam;
    const double rmplamalam = invarMass(plamcm, alamcm, plamalam);
    double ttplamcm, ttalamcm, ttp1alam, ttprecoil, ttplamalam;
    double ttmin, ttprimeplamcm, ttprimealamcm, ttprimep1alam,
        ttprimeprecoil, ttprimeplamalam;
    gett(pbeamcm, plamcm, ttplamcm, ttmin, ttprimeplamcm);
    gett(pbeamcm, alamcm, ttalamcm, ttmin, ttprimealamcm);
    gett(pbeamcm, ppalam, ttp1alam, ttmin, ttprimep1alam);
    gett(ptargetcm, p1cm, ttprecoil, ttmin, ttprimeprecoil);
    gett(pbeamcm, plamalam, ttplamalam, ttmin, ttprimeplamalam);

    FourVec ptemp;
    const FourVec& bachelor = (mech == 2) ? p2cm : p1cm;
    invarMass(plamalam, bachelor, ptemp);
    rmpalam = invarMass(bachelor, alamcm, ppalam);

    double rmiss, rmiss2;
    missingmass(psystem, ptemp, rmiss, rmiss2);
    im.hf1(9887, rmiss2, weight);
    const double rmissenergy = psystem.e - ptemp.e;
    im.hf1(9620, rmissenergy, weight);
    im.hf2(9621, rmiss2, rmissenergy, weight);

    double pmomplamcm, thplamcm, costhplamcm, phplamcm;
    double pmomalamcm, thalamcm, costhalamcm, phalamcm;
    double pmomplamlab, thplamlab, costhplamlab, phplamlab;
    double pmomalamlab, thalamlab, costhalamlab, phalamlab;
    double pmomplamalamcm, thplamalamcm, costhplamalamcm, phplamalamcm;
    getangles(plamcm, pmomplamcm, thplamcm, costhplamcm, phplamcm);
    getangles(alamcm, pmomalamcm, thalamcm, costhalamcm, phalamcm);
    getangles(plamlab, pmomplamlab, thplamlab, costhplamlab, phplamlab);
    getangles(alamlab, pmomalamlab, thalamlab, costhalamlab, phalamlab);
    getangles(plamalam, pmomplamalamcm, thplamalamcm, costhplamalamcm,
              phplamalamcm);

    im.hf1(9801, rmplam, weight);
    im.hf1(9802, rmalam, weight);
    im.hf1(9851, rmplamalam, weight);
    im.hf1(9852, rmpalam, weight);
    im.hf1(9860, rmpplam, weight);
    im.hf1(9856, -ttprimeplamalam, weight);
    im.hf1(9857, -ttprimep1alam, weight);
    im.hf1(9859, -ttprimeprecoil, weight);
    im.hf1(9861, -ttprimeplamcm, weight);
    im.hf1(9862, -ttprimealamcm, weight);
    im.hf2(9850, rmpplam, rmpalam, weight);
    im.hf2(9855, rmplamalam, rmpalam, weight);
    im.hf2(9854, rmpplam * rmpplam, rmpalam * rmpalam, weight);
    im.hf2(9849, rmplamalam * rmplamalam, rmpalam * rmpalam, weight);
    im.hf2(9625, rmplamalam, -ttprimeplamcm, weight);
    im.hf2(9626, rmplamalam, -ttprimealamcm, weight);
    im.hf2(9627, rmplamalam, -ttprimeprecoil, weight);
    im.hf2(9628, rmplamalam, -ttprimeplamalam, weight);
    im.hf2(9630, -ttplamalam, -ttprecoil, weight);
    im.hf2(9631, -ttprimeplamalam, -ttprimeprecoil, weight);

    im.hf1(28, ppbeam, weight);
    const Vec3 avertex = ev[k].rCreate;
    const Vec3 plampath{dvertex[0] - rvertex[0], dvertex[1] - rvertex[1],
                        dvertex[2] - rvertex[2]};
    const Vec3 alampath{avertex[0] - rvertex[0], avertex[1] - rvertex[1],
                        avertex[2] - rvertex[2]};
    const Vec3 dlampath{avertex[0] - dvertex[0], avertex[1] - dvertex[1],
                        avertex[2] - dvertex[2]};
    auto len3 = [](const Vec3& v) {
        return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    };
    const double plamlength = len3(plampath);
    const double alamlength = len3(alampath);
    const double alamlamlen = len3(dlampath);
    const double tauplam = plamlength * kMassLambda / (kClight * pmomplamlab);
    const double taualam = alamlength * kMassLambda / (kClight * pmomalamlab);
    im.hf1(604, rvertex[0], weight);
    im.hf1(605, rvertex[1], weight);
    im.hf1(606, rvertex[2], weight);
    im.hf1(607, dvertex[0], weight);
    im.hf1(608, dvertex[1], weight);
    im.hf1(609, dvertex[2], weight);
    im.hf1(610, avertex[0], weight);
    im.hf1(611, avertex[1], weight);
    im.hf1(612, avertex[2], weight);
    im.hf1(613, plamlength, weight);
    im.hf1(614, alamlength, weight);
    im.hf1(615, alamlamlen, weight);
    im.hf1(616, tauplam, weight);
    im.hf1(617, taualam, weight);
    im.hf2(618, taualam, tauplam, weight);

    im.hf2(9881, pmomplamcm, costhplamcm, weight);
    im.hf2(9883, pmomalamcm, costhalamcm, weight);
    im.hf2(9884, pmomplamalamcm, costhplamalamcm, weight);
    im.hf2(9869, costhplamcm, costhalamcm, weight);
    im.hf2(9880, pmomplamcm, pmomalamcm, weight);

    // Production-plane normal and the Lambda anti-Lambda decay normal.
    FourVec nhatdircm, pollamlam;
    crossprod(pbeamcm, plamalam, nhatdircm);
    crossprod(plamcm, alamcm, pollamlam);
    const double costh_planes =
        dot3(nhatdircm, pollamlam) / (nhatdircm.p() * pollamlam.p());
    const double planeangle = std::acos(costh_planes) * kRade;
    im.hf1(9885, costhplamcm, weight);
    im.hf1(9886, costhalamcm, weight);
    im.hf1(9888, costh_planes, weight);
    im.hf1(9889, planeangle, weight);

    // Pair rest frame.
    const Vec3 betapair = getbeta(plamalam);
    FourVec plampair = boost(plamcm, betapair);
    FourVec alampair = boost(alamcm, betapair);
    double pmomplampair, thplampair, costhplampair, phplpair;
    double pmomalampair, thalampair, costhalampair, phalpair;
    getangles(plampair, pmomplampair, thplampair, costhplampair, phplpair);
    getangles(alampair, pmomalampair, thalampair, costhalampair, phalpair);
    double costhplamdecay, thetaplamdecay, costhalamdecay, thetaalamdecay;
    getdecayangle(plampair, betapair, costhplamdecay, thetaplamdecay);
    getdecayangle(alampair, betapair, costhalamdecay, thetaalamdecay);
    im.hf1(9603, costhplamdecay, weight);
    im.hf1(9604, costhalamdecay, weight);
    im.hf1(9611, pmomplampair, weight);
    im.hf1(9612, pmomalampair, weight);
    const double pscale = 0.5;
    const double reweight = weight / pmomplampair * pscale;
    const double reweight2 = reweight / pmomplamalamcm;
    im.hf1(9616, pmomplampair, reweight);
    im.hf1(9617, pmomplampair, reweight2);
    im.hf1(9853, rmplamalam, reweight);
    im.hf1(9858, rmplamalam, reweight2);
    if (costhplamcm > 0.0) {
        im.hf1(9613, costhplamdecay, weight);
        im.hf1(9614, costhalamdecay, weight);
        im.hf1(9615, pmomplampair, weight);
    }

    // Van Hove and alternates.
    const double q1 = plamcm.pz;
    const double q2 = alamcm.pz;
    const double q3 = bachelor.pz;
    const int isector = vanHoveSector(q1, q2, q3);
    const double omega =
        std::atan2(-std::sqrt(3.0) * q1, q1 + 2.0 * q2) + M_PI;
    double pairmass = rmplamalam;
    if (isector == 1 || isector == 4) pairmass = rmpplam;
    else if (isector == 2 || isector == 5) pairmass = rmpalam;
    if (isector > 0) {
        const double qq = 0.81650 * std::sqrt(q1 * q1 + q2 * q2 + q3 * q3);
        const double cosv = std::cos(omega), sinv = std::sin(omega);
        im.hf2(9870, qq * cosv, qq * sinv, weight);
        im.hf2(9871, pairmass * cosv, pairmass * sinv, weight);
        const double rkinetic = ww - (2.0 * kMassLambda + kMassProton);
        if (rkinetic >= 0.0) {
            im.hf2(9872, rkinetic * cosv, rkinetic * sinv, weight);
        } else if (iloop == 1) {
            std::printf(" Negative available kinetic energy %g\n", rkinetic);
        }
        auto pt = [](const FourVec& p) {
            return std::sqrt(p.px * p.px + p.py * p.py);
        };
        const double ptrans = std::sqrt(pt(plamcm) * pt(plamcm) +
                                        pt(alamcm) * pt(alamcm) +
                                        pt(bachelor) * pt(bachelor));
        im.hf2(9873, ptrans * cosv, ptrans * sinv, weight);
    } else if (iloop == 1) {
        std::printf(" Cannot sort right... %g %g %g\n", q1, q2, q3);
    }
    im.hf2(9874, omega * kRade, pairmass, weight);

    // (Anti)proton in its parent hyperon's rest frame; direction cosines
    // against the {pollamlam, hyperon direction} axes.
    FourVec p1pair = boost(p1cm, betapair);
    FourVec p2pair = boost(p2cm, betapair);
    FourVec appair = boost(apcm, betapair);
    const Vec3 betatophyp = getbeta(plampair);
    const FourVec& lamProtonPair = (mech == 2) ? p1pair : p2pair;
    FourVec pphyp = boost(lamProtonPair, betatophyp);
    const Vec3 betatoahyp = getbeta(alampair);
    FourVec aphyp = boost(appair, betatoahyp);

    auto axes = [&](const Vec3& zdir, double* costh, const FourVec& probe) {
        FourVec zhat, yhat = pollamlam, xhat;
        zhat.px = zdir[0];
        zhat.py = zdir[1];
        zhat.pz = zdir[2];
        crossprod(yhat, zhat, xhat);
        const double xn = xhat.p(), yn = yhat.p(), zn = zhat.p();
        const double pm = probe.p();
        costh[0] = dot3(probe, xhat) / (pm * xn);
        costh[1] = dot3(probe, yhat) / (pm * yn);
        costh[2] = dot3(probe, zhat) / (pm * zn);
    };
    double cthp[3], ctha[3];
    axes(betatophyp, cthp, pphyp);
    axes(betatoahyp, ctha, aphyp);

    im.hf1(9601, cthp[1], weight);
    im.hf1(9602, ctha[1], weight);
    im.hf1(9605, cthp[0], weight);
    im.hf1(9608, ctha[0], weight);
    // NB: the Fortran fills both 9607 and 9610 with the ANTI-hyperon z
    // cosine (9606/9609 are booked but never filled); preserved as-is.
    im.hf1(9607, ctha[2], weight);
    im.hf1(9610, ctha[2], weight);

    if (iloop == 1) {
        auto& cs = im.corrsum;
        auto& cerr = im.correrr;
        for (int a = 0; a < 3; ++a) {
            cs[a][3] += cthp[a];
            cerr[a][3] += cthp[a] * cthp[a];
            cs[a][4] += ctha[a];
            cerr[a][4] += ctha[a] * ctha[a];
            for (int b = 0; b < 3; ++b) {
                cs[a][b] += ctha[a] * cthp[b];
                cerr[a][b] += (ctha[a] * cthp[b]) * (ctha[a] * cthp[b]);
            }
        }
    }
    for (int a = 0; a < 3; ++a) {
        im.hf2(9890, 4.0, a + 1.0, cthp[a]);
        im.hf2(9890, 5.0, a + 1.0, ctha[a]);
        for (int b = 0; b < 3; ++b)
            im.hf2(9890, a + 1.0, b + 1.0, ctha[a] * cthp[b]);
    }
    return true;  // first successful combination processed (goto 920)

                    }
                }
            }
        }
    }
    return true;
}

bool RootAnalyzer::ksksBlock(const Event& ev, double ppbeam, double weight,
                             int iloop, bool iprocess,
                             const std::vector<double>& acceptlist) {
    auto& im = *impl_;
    const Vec3 beta = getbeta(ev[0].p);
    FourVec pbeamlab{0.0, 0.0, ppbeam, 0.0, ppbeam};

    for (size_t i = 1; i < ev.size(); ++i) {
        if (ev[i].lund != 41) continue;
        for (size_t j = 1; j < ev.size(); ++j) {
            if (ev[j].lund != 17) continue;
            for (size_t k = j + 1; k < ev.size(); ++k) {
                if (ev[k].lund != 17) continue;
                for (size_t l = 1; l < ev.size(); ++l) {
                    if (ev[l].lund != -17) continue;
                    for (size_t m = l + 1; m < ev.size(); ++m) {
                        if (ev[m].lund != -17) continue;

    if (iloop == 1 && iprocess) ++im.ikeep;
    if (iloop == 3)
        weight = acceptlist[i] * acceptlist[j] * acceptlist[k] *
                 acceptlist[l] * acceptlist[m];

    double pmomlab, thplab, cthplab, phplab;
    double pip1momlab, th1p, cth1p, ph1p;
    double pip2momlab, th2p, cth2p, ph2p;
    double pim1momlab, th1m, cth1m, ph1m;
    double pim2momlab, th2m, cth2m, ph2m;
    getangles(ev[i].p, pmomlab, thplab, cthplab, phplab);
    getangles(ev[j].p, pip1momlab, th1p, cth1p, ph1p);
    getangles(ev[k].p, pip2momlab, th2p, cth2p, ph2p);
    getangles(ev[l].p, pim1momlab, th1m, cth1m, ph1m);
    getangles(ev[m].p, pim2momlab, th2m, cth2m, ph2m);
    im.hf1(10901, pmomlab, weight);
    im.hf1(10902, pip1momlab, weight);
    im.hf1(10903, pip2momlab, weight);
    im.hf1(10904, pim1momlab, weight);
    im.hf1(10905, pim2momlab, weight);
    im.hf1(10906, cthplab, weight);
    im.hf1(10907, cth1p, weight);
    im.hf1(10908, cth2p, weight);
    im.hf1(10909, cth1m, weight);
    im.hf1(10910, cth2m, weight);

    const bool gen = iloop == 1;
    auto pick = [&](size_t idx) { return gen ? ev[idx].p : ev[idx].pDecay; };
    FourVec protcm = boost(pick(i), beta);
    FourVec pip1cm = boost(pick(j), beta);
    FourVec pip2cm = boost(pick(k), beta);
    FourVec pim1cm = boost(pick(l), beta);
    FourVec pim2cm = boost(pick(m), beta);
    FourVec pbeamcm = boost(pbeamlab, beta);
    FourVec psystem = boost(ev[0].p, beta);

    double pmomcm, thpcm, cthpcm, phpcm;
    double pip1momcm, thpip1cm, cthpip1cm, phpip1cm;
    double pip2momcm, thpip2cm, cthpip2cm, phpip2cm;
    double pim1momcm, thpim1cm, cthpim1cm, phpim1cm;
    double pim2momcm, thpim2cm, cthpim2cm, phpim2cm;
    getangles(protcm, pmomcm, thpcm, cthpcm, phpcm);
    getangles(pip1cm, pip1momcm, thpip1cm, cthpip1cm, phpip1cm);
    getangles(pip2cm, pip2momcm, thpip2cm, cthpip2cm, phpip2cm);
    getangles(pim1cm, pim1momcm, thpim1cm, cthpim1cm, phpim1cm);
    getangles(pim2cm, pim2momcm, thpim2cm, cthpim2cm, phpim2cm);
    im.hf1(10911, pmomcm, weight);
    im.hf1(10912, pip1momcm, weight);
    im.hf1(10913, pip2momcm, weight);
    im.hf1(10914, pim1momcm, weight);
    im.hf1(10915, pim2momcm, weight);
    im.hf1(10916, cthpcm, weight);
    im.hf1(10917, cthpip1cm, weight);
    im.hf1(10918, cthpip2cm, weight);
    im.hf1(10919, cthpim1cm, weight);
    im.hf1(10920, cthpim2cm, weight);

    FourVec pks1cm, pks2cm;
    double rmks1 = invarMass(pip1cm, pim1cm, pks1cm);
    double rmks2 = invarMass(pip2cm, pim2cm, pks2cm);

    // Mimic the kinematic fit: reject bad kaon masses (iloop 2/3 only).
    if (iloop > 1) {
        if (rmks1 < 0.45 || rmks1 > 0.55) return false;
        if (rmks2 < 0.45 || rmks2 > 0.55) return false;
    }

    double pmomk1cm, thk1cm, costhk1cm, phk1cm;
    double pmomk2cm, thk2cm, costhk2cm, phk2cm;
    getangles(pks1cm, pmomk1cm, thk1cm, costhk1cm, phk1cm);
    getangles(pks2cm, pmomk2cm, thk2cm, costhk2cm, phk2cm);

    // Randomize (iloop 1) or sort by forwardness (iloop 2/3) the K0 labels.
    bool swap = false;
    if (iloop == 1) swap = im.rng.uniform() > 0.5;
    else swap = costhk1cm < costhk2cm;
    if (swap) {
        std::swap(pks1cm, pks2cm);
        std::swap(rmks1, rmks2);
    }
    getangles(pks1cm, pmomk1cm, thk1cm, costhk1cm, phk1cm);
    getangles(pks2cm, pmomk2cm, thk2cm, costhk2cm, phk2cm);

    FourVec pks12cm, ppks1cm, ppks2cm, ptemp;
    const double rmks12 = invarMass(pks1cm, pks2cm, pks12cm);
    const double rmpks1 = invarMass(pks1cm, protcm, ppks1cm);
    const double rmpks2 = invarMass(pks2cm, protcm, ppks2cm);
    invarMass(pks12cm, protcm, ptemp);
    double rmiss, rmiss2;
    missingmass(psystem, ptemp, rmiss, rmiss2);
    im.hf1(10887, rmiss2, weight);
    im.hf1(10888, rmiss2, weight);
    const double rmissenergy = psystem.e - ptemp.e;
    im.hf1(10620, rmissenergy, weight);
    im.hf2(10621, rmiss2, rmissenergy, weight);

    double pmomks12cm, thks12cm, costhks12cm, phks12cm;
    getangles(pks12cm, pmomks12cm, thks12cm, costhks12cm, phks12cm);
    double ttks1, tt1min, ttprimeks1, ttks2, tt2min, ttprimeks2;
    double ttks12, tt12min, ttprimeks12;
    gett(pbeamcm, pks1cm, ttks1, tt1min, ttprimeks1);
    gett(pbeamcm, pks2cm, ttks2, tt2min, ttprimeks2);
    gett(pbeamcm, pks12cm, ttks12, tt12min, ttprimeks12);

    im.hf1(10801, rmks1, weight);
    im.hf1(10802, rmks2, weight);
    im.hf1(10851, rmks12, weight);
    im.hf1(10854, rmks12, weight);
    im.hf1(10852, rmpks1, weight);
    im.hf1(10853, rmpks2, weight);
    im.hf1(10885, costhk1cm, weight);
    im.hf1(10886, costhk2cm, weight);
    im.hf2(10622, rmks1, rmks2, weight);
    im.hf2(10623, rmpks1, rmpks2, weight);
    im.hf2(10624, rmks12, rmpks2, weight);
    if (costhk2cm > 0.0) {
        im.hf1(10451, rmks12, weight);
        im.hf1(10452, rmpks1, weight);
    } else {
        im.hf1(10453, rmpks2, weight);
    }
    im.hf2(10889, costhk1cm, costhk2cm, weight);
    im.hf2(10921, rmpks2, -ttprimeks1, weight);
    im.hf2(10922, rmks12, -ttprimeks1, weight);
    im.hf2(10922, rmks12, -ttprimeks2, weight);
    im.hf2(10923, rmpks2, -ttprimeks12, weight);
    im.hf2(10924, rmks12, -ttprimeks12, weight);
    im.hf2(10931, rmpks2, costhk1cm, weight);
    im.hf2(10932, rmks12, costhk1cm, weight);
    im.hf2(10932, rmks12, costhk2cm, weight);
    im.hf2(10933, rmpks2, costhks12cm, weight);
    im.hf2(10934, rmks12, costhks12cm, weight);
    im.hf2(10941, -ttprimeks1, -ttprimeks2, weight);
    im.hf2(10942, costhk1cm, costhk2cm, weight);
    im.hf1(10951, costhk1cm, weight);
    im.hf1(10952, costhk2cm, weight);

    // Split the Dalitz-like plane into meson-like vs baryon-like regions.
    const double yy = 1.5 + (2.5 / 3.0) * (rmks12 - 1.0);
    if (rmpks2 > yy) im.hf1(10630, rmks12, weight);
    else im.hf1(10631, rmpks2, weight);

    im.hf1(28, ppbeam, weight);

    // Lab-frame kaons and vertexing.
    FourVec pks1lab, pks2lab;
    invarMass(pick(j), pick(l), pks1lab);
    invarMass(pick(k), pick(m), pks2lab);
    double pmomks1lab, thks1lab, cthks1lab, phks1lab;
    double pmomks2lab, thks2lab, cthks2lab, phks2lab;
    getangles(pks1lab, pmomks1lab, thks1lab, cthks1lab, phks1lab);
    getangles(pks2lab, pmomks2lab, thks2lab, cthks2lab, phks2lab);

    const Vec3 rvertex = ev[i].rCreate;
    const Vec3 dvertex = ev[j].rCreate;
    const Vec3 avertex = ev[k].rCreate;
    auto len3 = [](const Vec3& v) {
        return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    };
    const Vec3 v1{dvertex[0] - rvertex[0], dvertex[1] - rvertex[1],
                  dvertex[2] - rvertex[2]};
    const Vec3 v2{avertex[0] - rvertex[0], avertex[1] - rvertex[1],
                  avertex[2] - rvertex[2]};
    const Vec3 v12{avertex[0] - dvertex[0], avertex[1] - dvertex[1],
                   avertex[2] - dvertex[2]};
    const double pks1length = len3(v1);
    const double pks2length = len3(v2);
    const double pks1ks2len = len3(v12);
    const double tauks1 = pks1length * kMassK0 / (kClight * pmomks1lab);
    const double tauks2 = pks2length * kMassK0 / (kClight * pmomks2lab);
    im.hf1(1604, rvertex[0], weight);
    im.hf1(1605, rvertex[1], weight);
    im.hf1(1606, rvertex[2], weight);
    im.hf1(1607, dvertex[0], weight);
    im.hf1(1608, dvertex[1], weight);
    im.hf1(1609, dvertex[2], weight);
    im.hf1(1610, avertex[0], weight);
    im.hf1(1611, avertex[1], weight);
    im.hf1(1612, avertex[2], weight);
    im.hf1(1613, pks1length, weight);
    im.hf1(1614, pks2length, weight);
    im.hf1(1615, pks1ks2len, weight);
    im.hf1(1616, tauks1, weight);
    im.hf1(1617, tauks2, weight);
    im.hf2(1618, tauks1, tauks2, weight);
    im.hf1(1619, tauks1 + tauks2, weight);
    im.hf1(1622, pks1length + pks2length, weight);
    im.hf2(1623, pks1length + pks2length, tauks1 + tauks2, weight);

    // Van Hove.
    const double q1 = pks1cm.pz, q2 = pks2cm.pz, q3 = protcm.pz;
    const int isector = vanHoveSector(q1, q2, q3);
    const double omega =
        std::atan2(-std::sqrt(3.0) * q1, q1 + 2.0 * q2) + M_PI;
    if (isector > 0) {
        const double qq = 0.81650 * std::sqrt(q1 * q1 + q2 * q2 + q3 * q3);
        im.hf2(1620, qq * std::cos(omega), qq * std::sin(omega), weight);
        im.hf2(1621, rmpks1, omega * kRade, weight);
        im.hf2(1621, rmpks2, omega * kRade, weight);
        im.hf1(1624, omega * kRade, weight);
    } else if (iloop == 1) {
        std::printf(" Cannot sort right... %g %g %g\n", q1, q2, q3);
    }

                        }
                    }
                }
            }
    }
    return true;
}

void RootAnalyzer::finish() {
    auto& im = *impl_;
    const double alpha = 0.642, alphabar = -0.642;
    const double n = static_cast<double>(im.icount > 0 ? im.icount : 1);
    std::printf(" Analyzed %ld events; passed %ld events (%.2f%%)\n",
                im.icount, im.ikeep, 100.0 * im.ikeep / n);

    std::printf("\n Lambda anti-Lambda analysis: polarization results (%%):\n");
    auto& cs = im.corrsum;
    auto& ce = im.correrr;
    for (int ii = 0; ii < 3; ++ii) {
        for (int col : {3, 4}) {
            cs[ii][col] /= n;
            ce[ii][col] = std::sqrt(ce[ii][col] / n - cs[ii][col] * cs[ii][col]);
            cs[ii][col] *= (3.0 / alpha) * 100.0;
            ce[ii][col] *= (3.0 / alpha) / std::sqrt(n) * 100.0;
        }
    }
    std::printf("      Lambda Polarization (X,Y,Z): %8.2f%8.2f%8.2f   +-%8.2f%8.2f%8.2f\n",
                cs[0][3], cs[1][3], cs[2][3], ce[0][3], ce[1][3], ce[2][3]);
    std::printf(" anti-Lambda Polarization (X,Y,Z): %8.2f%8.2f%8.2f   +-%8.2f%8.2f%8.2f\n",
                cs[0][4], cs[1][4], cs[2][4], ce[0][4], ce[1][4], ce[2][4]);

    std::printf("\n Correlation Matrix:\n");
    const double rnorm = 9.0 / (alpha * alphabar);
    for (int ii = 0; ii < 3; ++ii) {
        for (int jj = 0; jj < 3; ++jj) {
            cs[ii][jj] /= n;
            ce[ii][jj] = std::sqrt(ce[ii][jj] / n - cs[ii][jj] * cs[ii][jj]);
            cs[ii][jj] *= rnorm * 100.0;
            ce[ii][jj] = std::abs(rnorm * ce[ii][jj] / std::sqrt(n) * 100.0);
        }
        std::printf("   %8.2f%8.2f%8.2f   +-%8.2f%8.2f%8.2f\n", cs[ii][0],
                    cs[ii][1], cs[ii][2], ce[ii][0], ce[ii][1], ce[ii][2]);
    }
    const double sf =
        25.0 * (1.0 + (cs[0][0] - cs[1][1] + cs[2][2]) / 100.0);
    std::printf("\n Singlet Fraction (%%) %8.2f\n\n", sf);

    std::printf(" Writing histogram file: %s\n", im.file->GetName());
    im.file->Write();
    im.file->Close();
}

}  // namespace mcgen
