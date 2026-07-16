#!/bin/sh
# Mech6 smoke validation WITH the ROOT analyzer enabled: generates
# gamma p -> p (p pbar) events, fills the mc_analyze histogram suite, and
# renders a review PDF of key spectra.
#
# Run from the repository root after building with ROOT:
#   sh validation/run_mech6_analyze.sh
# Evidence lands in validation/work/ (gitignored):
#   mech6a.root  - RAWEVENT/ACCEPTED/WEIGHTED histogram directories
#   mech6a.pdf   - multi-page plots of key histograms
set -e
BIN=${1:-build/src/app/mc_gen}
ROOTEXE=${ROOTEXE:-$HOME/CERN_ROOT/root_build/bin/root}
WORK=validation/work
mkdir -p "$WORK/genout"
cd "$WORK"

python3 - <<'EOF' > flux_test.ascii
for i in range(1000):
    e = 5.005 + 0.010 * i
    print(f"{e:.3f} {1000.0 if 6.3 < e < 11.7 else 0.0:.1f}")
EOF

# Enable the analyzer flag (line: "0 !Flag to select whether to call local
# event analyzer") in addition to the usual substitutions.
sed -e 's/TOTALSIZE/50000/g' -e 's/SAMPLESIZE/50000/g' \
    -e 's#GENDIR/#genout/#' \
    -e 's#/home/haoli/[^ ]*\.ascii#flux_test.ascii#' \
    -e 's/DYNCODE1/1.71/' -e 's/DYNCODE2/0.45/' -e 's/DYNCODE3/0.4/' \
    -e 's#^0[[:space:]]*!Flag to select#1	!Flag to select#' \
    ../gen_mech6_tCutOff.def > mech6a.def

../../"$BIN" mech6a.def > mc_gen_analyze.log 2>&1
cd ../..

"$ROOTEXE" -l -b -q "validation/plot_hists.C(\"$WORK/mech6a.root\", \"$WORK/mech6a.pdf\")"
echo "Evidence written: $WORK/mech6a.root and $WORK/mech6a.pdf"
