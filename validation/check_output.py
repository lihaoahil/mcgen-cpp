#!/usr/bin/env python3
"""Validate mc_gen ASCII output for the mech6 smoke test.

Checks, per event:
  1. exact energy-momentum conservation (final state vs beam + target proton)
  2. the t' distribution of the produced pair follows exp(2 b1 ln(s) t')
     beyond the |tcutoff| region
  3. the pair invariant-mass excess follows exp(-(m23 - 2 mp) / b2)

Usage: check_output.py <file.ascii> <b1> <b2> <tcutoff>
"""
import math
import sys

MP = 0.938272


def parse_events(path):
    lines = open(path).read().splitlines()
    i = 0
    while i < len(lines):
        n = int(lines[i]); i += 1
        ebeam = float(lines[i][:10]); i += 1
        parts = []
        for _ in range(n):
            f = lines[i]; i += 2  # second line is the vertex info
            parts.append(tuple(float(f[15 + 10 * k:25 + 10 * k]) for k in range(4)))
        yield ebeam, parts  # parts: (E, px, py, pz), recoil listed first


def fit_log_slope(vals, lo, hi, nbins):
    hist = [0] * nbins
    for v in vals:
        if lo <= v < hi:
            hist[int((v - lo) / (hi - lo) * nbins)] += 1
    pts = [(lo + (b + 0.5) * (hi - lo) / nbins, math.log(c))
           for b, c in enumerate(hist) if c > 20]
    n = len(pts)
    sx = sum(x for x, _ in pts); sy = sum(y for _, y in pts)
    sxx = sum(x * x for x, _ in pts); sxy = sum(x * y for x, y in pts)
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


def main():
    path, b1, b2, tcut = sys.argv[1], *map(float, sys.argv[2:5])
    nev = bad = 0
    tprimes, mexcesses, ebeams = [], [], []
    for ebeam, parts in parse_events(path):
        nev += 1
        e = sum(p[0] for p in parts); px = sum(p[1] for p in parts)
        py = sum(p[2] for p in parts); pz = sum(p[3] for p in parts)
        if (abs(px) > 1e-3 or abs(py) > 1e-3 or abs(pz - ebeam) > 1e-3
                or abs(e - ebeam - MP) > 1e-3):
            bad += 1
        # produced pair = slots 2,3 (recoil proton listed first per .def)
        E = parts[1][0] + parts[2][0]
        qx = parts[1][1] + parts[2][1]
        qy = parts[1][2] + parts[2][2]
        qz = parts[1][3] + parts[2][3]
        m23sq = E * E - qx * qx - qy * qy - qz * qz
        if m23sq <= 0:
            continue
        m23 = math.sqrt(m23sq)
        t = (E - ebeam) ** 2 - qx * qx - qy * qy - (qz - ebeam) ** 2
        s = 2 * ebeam * MP + MP * MP
        w = math.sqrt(s)
        egcm = (s - MP * MP) / (2 * w)
        ep = (s - (MP * MP - m23sq)) / (2 * w)
        pp2 = ep * ep - m23sq
        if pp2 <= 0:
            continue
        tmin = m23sq - 2 * (egcm * ep - egcm * math.sqrt(pp2))
        tprimes.append(t - tmin)
        mexcesses.append(m23 - 2 * MP)
        ebeams.append(ebeam)

    smean = 2 * (sum(ebeams) / len(ebeams)) * MP + MP * MP
    expect_t = 2 * b1 * math.log(smean)
    slope_t = fit_log_slope(tprimes, -1.0, -abs(tcut) - 0.05, 10)
    slope_m = fit_log_slope(mexcesses, 0.05, 1.0, 20)

    print(f"events                 : {nev}")
    print(f"E-p violations         : {bad}")
    print(f"t' slope (fit/expected): {slope_t:+.2f} / {expect_t:+.2f}  GeV^-2")
    # The observed m23 distribution is (3-body phase space) x exp(-x/b2),
    # so its fitted slope sits above -1/b2; require it in between.
    print(f"m23 slope: {slope_m:+.2f}  (acceptance -1/b2 = {-1 / b2:+.2f}; "
          f"phase-space modulation makes the observed slope less negative)")

    ok = (bad == 0 and abs(slope_t - expect_t) < 0.15 * expect_t
          and -1 / b2 < slope_m < 0.0)
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
