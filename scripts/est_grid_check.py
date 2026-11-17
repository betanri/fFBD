#!/usr/bin/env python3
"""Estimated shift times against exact numerical integration (validation test 1c).

The tree and the rates are held fixed (rates from the input's init_lambda/init_mu/init_psi, one value per interval), so the only free
parameters are the shift times. Their exact posterior and the exact marginal likelihood
    Z = integral of prior(times) x exp(lnFBD(tree | rates, times)) d(times)
are computed on a fine grid of ffbd's own FBD density (test option --eval-shifts: lnFBD of the start tree and rates for each line of
shift times; that density equals RevBayes dnBDSTP, test 1b). ffbd then samples the times (every other move switched off) and estimates Z
by stepping-stone with --ss-power fbd.

  est_grid_check.py one FFBD INPUT OUTDIR LO HI FIXED     one time in the window LO:HI, a second shift fixed at FIXED (> HI); 0.005 Myr grid
  est_grid_check.py two FFBD INPUT OUTDIR W               free:2 on [0, W] (two ordered times); 0.01 Myr grid (about 5.7 million cells)

The FBD density jumps where a shift time crosses a node or fossil age (the event moves to the other interval), so the grid must be fine:
on a 0.1 Myr grid the two-time posterior seemed off (KS p < 0.001); on the 0.01 Myr grid it agrees.
RevBayes 1.3.2 cannot be used for this check: during MCMC it does not recompute dnBDSTP when a stochastic timeline changes.
Standard library only.
"""
import sys, os, math, csv, subprocess, threading, pathlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sbc import ks_uniform

OFF = "node_slide,fossil_age,regraft_local,regraft_global,collapse_expand,root_slide,origin_slide,lambda_scale,mu_scale,psi_scale,lambda_mu_updown"

def evaluate(ffbd, inp, o, spec, extra, rows):
    """lnFBD for each tuple of shift times in rows (streamed through ffbd --eval-shifts)."""
    p = subprocess.Popen([ffbd, "run", inp, "--prefix", f"{o}/ev", "--mode", "control", "--rate-shifts", spec, "--eval-shifts", "-"] + extra,
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    def feed():
        for r in rows: p.stdin.write(" ".join(f"{t:.12g}" for t in r) + "\n")
        p.stdin.close()
    th = threading.Thread(target=feed); th.start()
    vals = [float(x) for x in p.stdout]; th.join(); p.wait()
    return vals

def run_ffbd(ffbd, inp, o, spec, extra, seeds=(1, 2)):
    logs, mls = [], []
    for s in seeds:
        subprocess.run([ffbd, "run", inp, "--prefix", f"{o}/mc{s}", "--mode", "control", "--seed", str(s), "--gens", "20000000", "--print", "2000",
                        "--tree-print", "20000000", "--burnin", "100000", "--no-moves", OFF, "--rate-shifts", spec] + extra, check=True, capture_output=True)
        logs.append(f"{o}/mc{s}.log")
        subprocess.run([ffbd, "run", inp, "--prefix", f"{o}/ss{s}", "--mode", "control", "--seed", str(s), "--burnin", "100000", "--no-moves", OFF,
                        "--rate-shifts", spec, "--ss", "50", "--ss-power", "fbd", "--ss-gens", "400000", "--ss-burn", "40000", "--ss-every", "40"] + extra,
                       check=True, capture_output=True)
        r = list(csv.DictReader(open(f"{o}/ss{s}.ml.tsv"), delimiter="\t")); mls.append(float(r[0]["stepping_stone"]))
    return logs, mls

def samples(logs, col):
    out = []
    for p in logs: out += [float(x[col]) for x in list(csv.DictReader(open(p), delimiter="\t"))[1:]]
    return out

def compare(name, lo, h, mass, smp):
    """mass[i]: posterior mass of the cell [lo + i h, lo + (i+1) h); smp: ffbd samples (thinned x5 for the KS test)."""
    S = sum(mass); cdf = []; c = 0.0
    for m in mass: c += m / S; cdf.append(c)
    n = len(mass); us = [cdf[min(max(int((x - lo) / h), 0), n - 1)] for x in smp[::5]]
    d, p = ks_uniform(us)
    mean = sum((lo + (i + 0.5) * h) * m for i, m in enumerate(mass)) / S
    print(f"{name}: exact mean {mean:.4f}, ffbd mean {sum(smp) / len(smp):.4f}; KS against the exact CDF D {d:.4f} p {p:.3f} (n {len(us)})")

if __name__ == "__main__":
    mode, ffbd, inp, outd = sys.argv[1:5]
    o = pathlib.Path(outd); o.mkdir(parents=True, exist_ok=True)
    if mode == "one":
        lo, hi, fx = map(float, sys.argv[5:8]); n = int(round((hi - lo) / 0.005)); h = (hi - lo) / n
        v = evaluate(ffbd, inp, str(o), f"{lo}:{hi},{fx}", [], ((lo + (i + 0.5) * h, fx) for i in range(n)))
        mx = max(v); w = [math.exp(x - mx) for x in v]
        lnZ = mx + math.log(sum(w) / n)                      # prior 1/(hi - lo): Z = mean of exp(lnFBD) over the window
        logs, mls = run_ffbd(ffbd, inp, str(o), f"{lo}:{hi},{fx}", [])
        compare("T1", lo, h, w, samples(logs, "shift_time[1]"))
        print(f"ln Z exact {lnZ:.4f} ({n} cells of {h:.4f} Myr); ffbd stepping-stone " + " ".join(f"{x:.4f}" for x in mls))
    elif mode == "two":
        W = float(sys.argv[5]); n = int(round(W / 0.01)); h = W / n; ex = ["--shift-window", f"0,{W}"]
        pairs = ((i, j) for i in range(n) for j in range(i + 1, n))
        v = evaluate(ffbd, inp, str(o), "free:2", ex, (((i + 0.5) * h, (j + 0.5) * h) for i in range(n) for j in range(i + 1, n)))
        mx = max(v); m1 = [0.0] * n; m2 = [0.0] * n; tot = 0.0
        for (i, j), x in zip(pairs, v):
            e = math.exp(x - mx); m1[i] += e; m2[j] += e; tot += e
        lnZ = mx + math.log(tot * h * h * 2 / W ** 2)        # prior 2/W^2 on t1 < t2
        logs, mls = run_ffbd(ffbd, inp, str(o), "free:2", ex)
        compare("T1", 0.0, h, m1, samples(logs, "shift_time[1]")); compare("T2", 0.0, h, m2, samples(logs, "shift_time[2]"))
        print(f"ln Z grid {lnZ:.4f} ({len(v)} cells of {h:.4f} Myr; cells on the diagonal, {n * h * h / 2 / (W * W / 2):.3%} of the area, left out); "
              f"ffbd stepping-stone " + " ".join(f"{x:.4f}" for x in mls))
