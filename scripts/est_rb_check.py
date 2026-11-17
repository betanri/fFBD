#!/usr/bin/env python3
"""Skyline FBD on a fixed tree: ffbd against RevBayes 1.3.2 (posterior of the rates and marginal likelihood). Validation test 1c.

The tree (topology, node ages, fossil ages, sampled ancestors) is held fixed in both programs, so the free parameters are the per-interval
rates lambda_i, mu_i, psi_i ~ Exp(10). The target is p(rates | tree) proportional to prior(rates) dnBDSTP(tree | rates, timeline), with
RevBayes dnBDSTP(rootAge|originAge, ..., timeline=...) clamped to the tree. The marginal likelihood p(tree) is estimated by ffbd
--ss-power fbd (tree moves switched off) and by RevBayes powerPosterior + steppingStoneSampler.

Shift times must be fixed here (--shifts T1,T2,...). RevBayes 1.3.2 cannot check estimated shift times: during MCMC it does not
recompute dnBDSTP when a stochastic timeline changes (the likelihood stays at its start value), so the times sample their prior.
Estimated times are checked against exact numerical integration instead (est_grid_check.py).

  est_rb_check.py gen INPUT OUTDIR --ffbd FFBD --shifts T1,T2 [--gens G] [--thin T] [--ss-stones K] [--ss-gens G]
                      [--rb-gens G] [--rb-ss-gens G] [--seeds 2] [--redo-rb]
      writes OUTDIR/run.sh (ffbd MCMC and SS for each seed, then RevBayes; completed steps are skipped) and OUTDIR/rb_S.Rev
  est_rb_check.py compare OUTDIR
      per rate: posterior mean, sd and 5/50/95% quantiles in both programs and the difference of means in units of its Monte Carlo
      standard error (batch means); the marginal likelihoods of every seed.
"""
import sys, os, pathlib, re, math, subprocess, statistics as st
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sky_density_check import read_dump, read_input, tip_ages

TREE_MOVES = "node_slide,fossil_age,regraft_local,regraft_global,collapse_expand,root_slide,origin_slide"

REV = r'''setOption("outputPrecision", "15", FALSE)
seed(@SEED@)
taxa <- readTaxonData("@DIR@/taxa.tsv")
tr <- readTrees("@DIR@/tree.nwk", treetype="clock")[1]
for (i in 1:@NI@) {
  lam[i] ~ dnExponential(10)
  mu[i] ~ dnExponential(10)
  psi[i] ~ dnExponential(10)
}
mi = 0
for (i in 1:@NI@) {
  moves[++mi] = mvScale(lam[i], lambda=0.5, weight=3)
  moves[++mi] = mvScale(mu[i], lambda=0.5, weight=3)
  moves[++mi] = mvScale(psi[i], lambda=0.5, weight=3)
  moves[++mi] = mvUpDownScale(lambda=0.3, weight=2)
  moves[mi].addVariable(lam[i], up=TRUE)
  moves[mi].addVariable(mu[i], up=TRUE)
}
@TIMES@
fbd ~ dnBDSTP(@COND@=@CONDV@, lambda=lam, mu=mu, psi=psi, rho=@RHO@, r=0, timeline=tl, taxa=taxa)
fbd.clamp(tr)
mymodel = model(lam)
monitors[1] = mnModel(filename="@DIR@/rb_@SEED@.log", printgen=@THIN@)
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.burnin(generations=@BURN@, tuningInterval=200)
mymcmc.run(generations=@GENS@)
ppmon[1] = mnModel(filename="@DIR@/rbpp_@SEED@.samples", printgen=@SSTHIN@)   # own monitor: keep the MCMC log free of power-posterior samples
pp = powerPosterior(mymodel, moves, ppmon, "@DIR@/rb_@SEED@.pp", cats=@CATS@, alpha=0.3, sampleFreq=@SSTHIN@)
pp.burnin(generations=@BURN@, tuningInterval=200)
pp.run(generations=@SSGENS@)
ss = steppingStoneSampler(file="@DIR@/rb_@SEED@.pp", powerColumnName="power", likelihoodColumnName="likelihood")
print("RBSS", ss.marginal(), sep="\t")
q()
'''

def times_block(spec, window):
    """Rev code for the shift times (timeline tl) and their moves; returns (code, K, names of the logged times)."""
    if spec.startswith("free:"):
        K = int(spec[5:]); lo, hi = window
        code = [f"for (k in 1:{K}) {{ t[k] ~ dnUniform({lo}, {hi}); t[k].setValue({lo} + k * ({hi} - {lo}) / {K + 1}) }}",
                f"for (k in 1:{K}) {{ moves[++mi] = mvSlide(t[k], delta={(hi - lo) / 10}, weight=5) }}",
                "tl := sort(t)"]
        return "\n".join(code), K
    items = spec.split(",")
    parts, code = [], []
    for j, it in enumerate(items):
        if ":" in it:
            lo, hi = map(float, it.split(":"))
            code.append(f"t{j} ~ dnUniform({lo}, {hi}); t{j}.setValue({(lo + hi) / 2})")
            code.append(f"moves[++mi] = mvSlide(t{j}, delta={(hi - lo) / 10}, weight=5)")
            parts.append(f"t{j}")
        else:
            parts.append(it)
    code.append("tl := v(" + ",".join(parts) + ")")
    return "\n".join(code), len(items)

def gen(a):
    o = pathlib.Path(a.outdir).resolve(); o.mkdir(parents=True, exist_ok=True)
    if a.redo_rb:
        for f in list(o.glob("rb_*")) + list(o.glob("rbpp_*")): f.unlink()
    rho, cond, types = read_input(a.input)
    win = [float(x) for x in a.window.split(",")] if a.window else None
    # the fixed tree: one ffbd proposal with every tree move switched off, dumped
    cmd = [a.ffbd, "run", a.input, "--mode", "control", "--prefix", str(o / "tree0"), "--seed", "1", "--gens", "1", "--print", "1",
           "--dump", str(o / "tree0.dump"), "--dump-every", "1", "--no-moves", TREE_MOVES]
    subprocess.run(cmd, check=True, capture_output=True)
    r = read_dump(o / "tree0.dump")[0]
    d = tip_ages(r["nwk"]); ext = [d[t] for t in d if types[t] == "E"]; h = max(ext)
    with open(o / "taxa.tsv", "w") as f:
        f.write("taxon\tmin_age\tmax_age\n")
        for t in d:
            age = max(0.0, h - d[t])
            if types[t] == "E": f.write(f"{t}\t0.0\t0.0\n")
            else: f.write(f"{t}\t{max(0.0, age - 1e-7):.15g}\t{age + 1e-7:.15g}\n")
    (o / "tree.nwk").write_text(r["nwk"] + "\n")
    if a.shifts.startswith("free:") or ":" in a.shifts:
        sys.exit("RevBayes 1.3.2 does not recompute dnBDSTP when a stochastic timeline moves; use fixed --shifts here and est_grid_check.py for estimated times")
    tcode, K = times_block(a.shifts, win)
    nI = K + 1
    sh = ["#!/bin/sh", "set -e", f"cd {o}"]
    ff = [a.ffbd, "run", os.path.abspath(a.input), "--mode", "control", "--no-moves", TREE_MOVES, "--rate-shifts", a.shifts]
    if win: ff += ["--shift-window", a.window]
    for s in range(1, a.seeds + 1):
        seed = 100 + s
        sh.append(" ".join(ff + ["--prefix", f"ff_{seed}", "--seed", str(seed), "--gens", str(a.gens), "--print", str(a.thin),
                                 "--tree-print", str(a.gens), "--burnin", str(a.gens // 10)]) + f" > ff_{seed}.out")
        sh[-1] = f"grep -q 'RUN done' ff_{seed}.out 2>/dev/null || " + sh[-1]
        sh.append(" ".join(ff + ["--prefix", f"ffss_{seed}", "--seed", str(seed), "--burnin", str(a.ss_gens), "--ss", str(a.ss_stones),
                                 "--ss-power", "fbd", "--ss-gens", str(a.ss_gens), "--ss-burn", str(a.ss_gens // 5),
                                 "--ss-every", str(max(1, a.ss_gens // 2000))]) + f" > ffss_{seed}.out")
        sh[-1] = f"[ -s ffss_{seed}.ml.tsv ] || " + sh[-1]
        txt = (REV.replace("@DIR@", str(o)).replace("@SEED@", str(seed)).replace("@NI@", str(nI)).replace("@TIMES@", tcode)
               .replace("@COND@", "originAge" if cond == "origin" else "rootAge")
               .replace("@CONDV@", r["cond"] if cond == "origin" else "tr.rootAge()").replace("@RHO@", rho)
               .replace("@THIN@", str(a.rb_thin)).replace("@BURN@", str(a.rb_gens // 10)).replace("@GENS@", str(a.rb_gens))
               .replace("@CATS@", str(a.ss_stones)).replace("@SSGENS@", str(a.rb_ss_gens)).replace("@SSTHIN@", str(max(1, a.rb_ss_gens // 1000))))
        (o / f"rb_{seed}.Rev").write_text(txt)
        sh.append(f"grep -q RBSS rb_{seed}.out 2>/dev/null || rb rb_{seed}.Rev > rb_{seed}.out 2>&1")
    (o / "run.sh").write_text("\n".join(sh) + "\n")
    print(f"wrote {o}/run.sh: K = {K} shifts ({a.shifts}), tree N = {len(d)}, cond {cond}")

def read_tsv(p, skip_first=False):
    with open(p) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in f if l.strip()]
    if skip_first: rows = rows[1:]
    return hdr, rows

def mcse(x, nb=20):
    n = len(x) // nb
    if n < 2: return float("nan")
    m = [st.mean(x[i * n:(i + 1) * n]) for i in range(nb)]
    return st.stdev(m) / math.sqrt(nb)

def q(x, p):
    s = sorted(x); return s[min(len(s) - 1, int(p * len(s)))]

def compare(outd):
    o = pathlib.Path(outd)
    ffl = sorted(o.glob("ff_*.log")); rbl = sorted(o.glob("rb_*.log"))
    if not ffl or not rbl: sys.exit("no logs")
    # ffbd columns: speciation_rate[i], extinction_rate[i], psi[i], shift_time[k]; RevBayes: lam[i], mu[i], psi[i], t[k] or tK
    F, R = {}, {}
    for p in ffl:
        h, rows = read_tsv(p, skip_first=True)
        for j, c in enumerate(h):
            mm = re.match(r"(speciation_rate|extinction_rate|psi|shift_time)\[(\d+)\]$", c)
            if mm:
                nm = {"speciation_rate": "lambda", "extinction_rate": "mu", "psi": "psi", "shift_time": "T"}[mm.group(1)] + mm.group(2)
                F.setdefault(nm, []).extend(float(r[j]) for r in rows)
    for p in rbl:
        h, rows = read_tsv(p)
        rows = rows[len(rows) // 10:]
        tcols = [j for j, c in enumerate(h) if re.match(r"t(\[\d+\]|\d+)$", c)]
        for j, c in enumerate(h):
            mm = re.match(r"(lam|mu|psi)\[(\d+)\]$", c)
            if mm: R.setdefault({"lam": "lambda", "mu": "mu", "psi": "psi"}[mm.group(1)] + mm.group(2), []).extend(float(r[j]) for r in rows)
        if tcols:
            ts = [sorted(float(r[j]) for j in tcols) for r in rows]
            fixedT = []
            for k in range(len(tcols)): R.setdefault(f"T{k + 1}", []).extend(t[k] for t in ts)
    # with fixed and windowed times mixed, RevBayes logs only the windowed ones: align by ffbd's estimated columns
    print("param\tffbd_mean\trb_mean\tz_mean\tffbd_sd\trb_sd\tffbd_q05/q50/q95\trb_q05/q50/q95")
    worst = 0
    est_ff = [k for k in F if k.startswith("T") and len(set(round(v, 12) for v in F[k][:200])) > 1]
    tmap = {}
    if [k for k in R if k.startswith("T")] and len(est_ff) == len([k for k in R if k.startswith("T")]):
        tmap = dict(zip(sorted([k for k in R if k.startswith("T")]), sorted(est_ff)))
    for k in sorted(F):
        rk = k
        if k.startswith("T"):
            inv = {v: kk for kk, v in tmap.items()}
            if k not in inv: continue
            rk = inv[k]
        if rk not in R: continue
        a, b = F[k], R[rk]
        se = math.sqrt(mcse(a) ** 2 + mcse(b) ** 2); z = (st.mean(a) - st.mean(b)) / se if se > 0 else float("nan")
        worst = max(worst, abs(z))
        print(f"{k}\t{st.mean(a):.5g}\t{st.mean(b):.5g}\t{z:+.2f}\t{st.stdev(a):.4g}\t{st.stdev(b):.4g}\t"
              f"{q(a, .05):.4g}/{q(a, .5):.4g}/{q(a, .95):.4g}\t{q(b, .05):.4g}/{q(b, .5):.4g}/{q(b, .95):.4g}")
    print(f"largest |z| of the posterior means: {worst:.2f} ({len(F)} parameters; |z| < 3 expected)")
    ml_ff = []
    for p in sorted(o.glob("ffss_*.ml.tsv")):
        h, rows = read_tsv(p); ml_ff.append(float(rows[0][h.index("stepping_stone")]))
    ml_rb = []
    for p in sorted(o.glob("rb_*.out")):
        for l in open(p, errors="replace"):
            if l.startswith("RBSS"): ml_rb.append(float(l.split("\t")[1]))
    print("marginal likelihood (ln) ffbd stepping-stone:", " ".join(f"{x:.4f}" for x in ml_ff))
    print("marginal likelihood (ln) RevBayes stepping-stone:", " ".join(f"{x:.4f}" for x in ml_rb))
    if ml_ff and ml_rb: print(f"difference of means ffbd - RevBayes: {st.mean(ml_ff) - st.mean(ml_rb):+.4f}")

if __name__ == "__main__":
    if sys.argv[1] == "compare": compare(sys.argv[2]); sys.exit(0)
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd"); ap.add_argument("input"); ap.add_argument("outdir")
    ap.add_argument("--ffbd", required=True); ap.add_argument("--shifts", required=True); ap.add_argument("--window", default="")
    ap.add_argument("--gens", type=int, default=4000000); ap.add_argument("--thin", type=int, default=1000)
    ap.add_argument("--ss-stones", type=int, default=50); ap.add_argument("--ss-gens", type=int, default=200000)
    ap.add_argument("--rb-gens", type=int, default=200000); ap.add_argument("--rb-thin", type=int, default=50)
    ap.add_argument("--rb-ss-gens", type=int, default=20000); ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--redo-rb", action="store_true", help="remove RevBayes outputs first, so run.sh reruns only RevBayes")
    gen(ap.parse_args())
