#!/bin/sh
# Mech6 smoke validation: gamma p -> p (p pbar) with b1=1.71 GeV^-2,
# b2=0.45 GeV, tcutoff=+0.4 GeV^2, 50k events over a flat synthetic flux.
# Run from the repository root after building:  sh validation/run_mech6.sh
set -e
BIN=${1:-build/src/app/mc_gen}
WORK=validation/work
mkdir -p "$WORK/genout"
cd "$WORK"

# Synthetic flat beam flux, 10 MeV bins (generator only needs the shape).
python3 - <<'EOF' > flux_test.ascii
for i in range(1000):
    e = 5.005 + 0.010 * i
    print(f"{e:.3f} {1000.0 if 6.3 < e < 11.7 else 0.0:.1f}")
EOF

sed -e 's/TOTALSIZE/50000/g' -e 's/SAMPLESIZE/50000/g' \
    -e 's#GENDIR/#genout/#' \
    -e 's#/home/haoli/[^ ]*\.ascii#flux_test.ascii#' \
    -e 's/DYNCODE1/1.71/' -e 's/DYNCODE2/0.45/' -e 's/DYNCODE3/0.4/' \
    ../gen_mech6_tCutOff.def > mech6.def

# Run inside the work dir: the .def references flux/output paths relatively.
../../"$BIN" mech6.def > mc_gen.log 2>&1
cd ../..
python3 validation/check_output.py "$WORK/genout/mech6_0001.ascii" 1.71 0.45 0.4
