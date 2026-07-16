# mcgen-cpp

A standalone, dependency-light C++20 remake of **MC_GEN**, the Monte Carlo
event generator for sequential N-body decays in photoproduction reactions
(GlueX / Hall D), currently part of
[JeffersonLab/halld_sim](https://github.com/JeffersonLab/halld_sim/tree/master/src/programs/Simulation/MC_GEN). This work is carried out in the context of the
[GlueX Collaboration](https://gluex.org) (Jefferson Lab, Hall D), whose
simulation ecosystem hosts the original MC_GEN.

MC_GEN was originally developed and maintained in Fortran by **Reinhard Schumacher**
at Carnegie Mellon University over three decades (1991–2020). In later applications, **Hao Li** contributed (while at Carnegie Mellon) to the phenomenological modeling of t-channel single and double-Regge exchange processes using MC_GEN-based workflows, with further details available in his [Ph.D. thesis](https://www.osti.gov/biblio/2301675). It was later integrated into the GlueX simulation pipeline through [gluex_MCwrapper](https://github.com/JeffersonLab/gluex_MCwrapper). This repository contains a modern C++ remake by **Hao Li**, intended to preserve the
original physics logic while improving readability, maintainability, and
extensibility for future GlueX analysis workflows.


MC_GEN and related sequential-decay Monte Carlo tools have been used in a number
of CLAS and GlueX photoproduction analyses, including recent baryon--antibaryon
photoproduction studies as well as earlier hyperon and meson production
measurements. Representative publications include:

- H. Li, R. A. Schumacher, et al. (GlueX Collaboration), "Baryon--Antibaryon Photoproduction Cross Sections off the Proton," [*Phys. Rev. C* **113**, 045207 (2026)](https://doi.org/10.1103/h197-l9dw), [arXiv:2510.26890 [hep-ex]](https://arxiv.org/abs/2510.26890).

- D. H. Ho, R. A. Schumacher, et al. (CLAS Collaboration), "Beam-Target Helicity Asymmetry $E$ in $K^0\Lambda$ and $K^0\Sigma^0$ Photoproduction on the Neutron," [*Phys. Rev. C* **98**, 045205 (2018)](https://doi.org/10.1103/PhysRevC.98.045205), [arXiv:1805.04561 [nucl-ex]](https://arxiv.org/abs/1805.04561).

- D. Ho, P. Peng, et al. (CLAS Collaboration), "Beam-Target Helicity Asymmetry for $\vec{\gamma}\vec{n} \to \pi^-p$ in the $N^\ast$ Resonance Region," [*Phys. Rev. Lett.* **118**, 242002 (2017)](https://doi.org/10.1103/PhysRevLett.118.242002), [arXiv:1705.04713 [nucl-ex]](https://arxiv.org/abs/1705.04713).

- R. Dickson, R. A. Schumacher, et al. (CLAS Collaboration), "Photoproduction of the $f_1(1285)$ Meson," [*Phys. Rev. C* **93**, 065202 (2016)](https://doi.org/10.1103/PhysRevC.93.065202), [arXiv:1604.07425 [nucl-ex]](https://arxiv.org/abs/1604.07425).

- K. Moriya et al. (CLAS Collaboration), "Spin and Parity Measurement of the $\Lambda(1405)$ Baryon," [*Phys. Rev. Lett.* **112**, 082004 (2014)](https://doi.org/10.1103/PhysRevLett.112.082004), [arXiv:1402.2296 [hep-ex]](https://arxiv.org/abs/1402.2296).

- K. Moriya et al. (CLAS Collaboration), "Differential Photoproduction Cross Sections of the $\Sigma^0(1385)$, $\Lambda(1405)$, and $\Lambda(1520)$," [*Phys. Rev. C* **88**, 045201 (2013)](https://doi.org/10.1103/PhysRevC.88.045201); Publisher's Note: [*Phys. Rev. C* **88**, 049902 (2013)](https://doi.org/10.1103/PhysRevC.88.049902), [arXiv:1305.6776 [nucl-ex]](https://arxiv.org/abs/1305.6776).

- R. A. Schumacher and K. Moriya, "Isospin Decomposition of the Photoproduced $\Sigma\pi$ System Near the $\Lambda(1405)$," [*Nucl. Phys. A* **914**, 51 (2013)](https://doi.org/10.1016/j.nuclphysa.2013.03.003), [arXiv:1303.0860 [nucl-ex]](https://arxiv.org/abs/1303.0860).

- K. Moriya et al. (CLAS Collaboration), "Measurement of the $\Sigma\pi$ Photoproduction Line Shapes Near the $\Lambda(1405)$," [*Phys. Rev. C* **87**, 035206 (2013)](https://doi.org/10.1103/PhysRevC.87.035206), [arXiv:1301.5000 [nucl-ex]](https://arxiv.org/abs/1301.5000).

- J. W. C. McNabb, R. A. Schumacher, L. Todor, et al. (CLAS Collaboration), "Hyperon Photoproduction in the Nucleon Resonance Region," [*Phys. Rev. C* **69**, 042201 (2004)](https://doi.org/10.1103/PhysRevC.69.042201), [arXiv:nucl-ex/0305028](https://arxiv.org/abs/nucl-ex/0305028).


## Purpose

The goal of this project is to provide a clean, portable, and maintainable implementation of the origin mc_gen in modern C++:

- generate sequential N-body decay chains for photoproduction reactions;
- preserve the physics behavior of the legacy MC_GEN code;
- make the code easier to inspect, validate, and extend;
- support future GlueX analyses that may need lightweight standalone event
  generation tools.


## Build and run

```sh
cmake -B build            # add -DMCGEN_WITH_ROOT=OFF to skip the analysis module
cmake --build build
ctest --test-dir build
build/src/app/mc_gen my_config.def
```

Requires CMake ≥ 3.20 and a C++20 compiler. The core generator has no
external dependencies; the (in progress) analysis module uses ROOT.

## Layout

- `src/core` — generator: four-vectors/kinematics (`kinpack2`, `carrot`),
  N-body phase space (`kin3body*`, GENBOD), lineshape samplers, beam
  generation, production dynamics, `.def` parser, event loop, ASCII writer
- `src/analysis` — ROOT-based histogramming (port of `mc_analyze`; HBOOK
  IDs become `h<ID>` TH1D/TH2D in `RAWEVENT`/`ACCEPTED`/`WEIGHTED`
  directories of a `.root` file)
- `tests` — physics unit tests (momentum conservation, distribution shapes)
- `validation` — end-to-end smoke validation against known physics

## Validation example

`sh validation/run_mech6.sh` generates 50,000 events of
`gamma p -> p (p pbar)` via "Mechanism 6" (t-slope b1 = 1.71 GeV^-2,
pair-mass slope b2 = 0.45 GeV, t-cutoff 0.4 GeV^2) over a flat synthetic
flux, then checks the output file:

```
events                 : 50000
E-p violations         : 0
t' slope (fit/expected): +9.50 / +9.84  GeV^-2
m23 slope: -1.79  (acceptance -1/b2 = -2.22; phase-space modulation ...)
PASS
```

Energy–momentum is conserved exactly in every event, and the fitted t'
slope matches the configured `2 b1 ln(s)` within statistics.

`sh validation/run_mech6_analyze.sh` runs the same configuration with the
ROOT analyzer enabled and renders review evidence into the gitignored
`validation/work/` directory: `mech6a.root` (full histogram suite) and
`mech6a.pdf` (key spectra: beam, IM(p2 pbar), t', cos theta_cm, Dalitz,
van Hove, acceptance-selected copies).

## Intentional deviations from the Fortran

- `double` precision throughout (the Fortran used single-precision `REAL`)
- `std::mt19937_64` replaces the GEANT3/RANECU random-number machinery;
  results are statistically equivalent, not bit-identical
- The truncated-Cauchy sampler (`bwran`) uses the analytic inverse CDF
  instead of a 1001-bin lookup table
- Dynamics lines in `.def` files with only 4 values (pre-04/2020 format)
  are accepted with a warning and `dyncode3 = 0`, instead of the Fortran's
  list-directed read spilling into the next line
- Copy-pasted Fortran subroutines that existed only to duplicate static
  state (`relativisticbwran1/2`) are single C++ classes instantiated per use

## Status

Ported and validated: full generator (`mc_gen.F`, `mc_propagate.F`, all
kinematics modules) and the `mc_analyze` analyzer with ROOT histogram
output (p pbar p, Lambda anti-Lambda, and Ks Ks p analysis blocks, GlueX
acceptance, polarization-correlation summary). Planned: distribution-level
comparison against reference output from the original Fortran/CERNLIB
build.
