#!/usr/bin/env python3
"""Score the partition schemes of test 5f (codon-position rates, sim_codons.py) on the evaluation replicates 51-100 against the true node ages, next to the
original arms on the same replicates. Same scoring as sim_genes_summary.py (first 10% dropped; 95% interval, posterior median).
  codon_test_summary.py SIM2DIR SIM1DIR [FIRST LAST] [OUT.tsv]
"""
import csv, os, statistics as st, sys

D2, D1 = sys.argv[1], sys.argv[2]
FIRST, LAST = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (51, 100)
OUT = sys.argv[5] if len(sys.argv) > 5 else None
ARMS = [('exact', '{d1}/{r}/exact.log', 'model-consistent data (b ~ N(true, sqrt(b/2000)))'),
        ('cod_unpart', '{d2}/{r}/codon/cod_unpart.log', '23 codon genes, one GTR+G4 for all sites'),
        ('cod_gene', '{d2}/{r}/codon/cod_gene.log', '23 codon genes, partitioned by gene (23)'),
        ('cod_codon', '{d2}/{r}/codon/cod_codon.log', '23 codon genes, partitioned by codon position (3)'),
        ('cod_genecodon', '{d2}/{r}/codon/cod_genecodon.log', '23 codon genes, gene x codon position (69, true scheme)'),
        ('cod_mfp', '{d2}/{r}/codon/cod_mfp.log', '23 codon genes, gene x codon merged by ModelFinder (MFP+MERGE)')]
BINS = [(0, 5), (5, 10), (10, 20), (20, 40), (40, 60), (60, 101)]


def q(x, a):
    i = a * (len(x) - 1); lo = int(i)
    return x[lo] + (x[min(lo + 1, len(x) - 1)] - x[lo]) * (i - lo)


rows, reps = [], []
for k in range(FIRST, LAST + 1):
    r = f'r{k}'
    logs = {a: p.format(d1=D1, d2=D2, r=r) for a, p, _ in ARMS}
    if not all(os.path.exists(f) and sum(1 for _ in open(f)) > 1000 for f in logs.values()):
        continue
    reps.append(r)
    truth = {l.split('\t')[0]: float(l.split('\t')[1]) for l in open(f'{D2}/{r}/truth.tsv').read().splitlines()[1:]}
    for a, f in logs.items():
        lg = list(csv.DictReader(open(f), delimiter='\t'))
        lg = lg[len(lg) // 10:]
        for p, t in truth.items():
            if not p.startswith('age_n') or p not in lg[0]:
                continue
            x = sorted(float(z[p]) for z in lg)
            lo, med, hi = q(x, .025), q(x, .5), q(x, .975)
            rows.append({'rep': r, 'arm': a, 'node': p[4:], 'true': t, 'median': med, 'lo95': lo, 'hi95': hi,
                         'in95': int(lo <= t <= hi), 'relerr': (med - t) / t, 'abserr': abs(med - t), 'width': hi - lo})
if OUT:
    with open(OUT, 'w') as fh:
        w = csv.DictWriter(fh, rows[0].keys(), delimiter='\t'); w.writeheader(); w.writerows(rows)
print(f'evaluation replicates {FIRST}-{LAST} with all arms: {len(reps)}; nodes per arm: {sum(1 for x in rows if x["arm"] == "exact")}')
print(f'{"arm":15} {"cov95":>6} {"median rel err":>14} {"median |err| Myr":>16} {"median width Myr":>16}')
for a, _, desc in ARMS:
    s = [x for x in rows if x['arm'] == a]
    print(f'{a:15} {st.fmean(x["in95"] for x in s):6.3f} {100 * st.median(x["relerr"] for x in s):13.1f}% '
          f'{st.median(x["abserr"] for x in s):16.2f} {st.median(x["width"] for x in s):16.2f}   {desc}')
print('\nby true node age (Ma): coverage / median relative error / median absolute error (Myr) / median width (Myr)')
print(f'{"arm":15} ' + ' '.join(f'{f"{lo}-{hi if hi < 100 else 100}":>24}' for lo, hi in BINS))
for a, _, _ in ARMS:
    out = []
    for lo, hi in BINS:
        s = [x for x in rows if x['arm'] == a and lo <= x['true'] < hi]
        out.append(f'{st.fmean(x["in95"] for x in s):.3f}/{100 * st.median(x["relerr"] for x in s):+.0f}%/{st.median(x["abserr"] for x in s):.2f}/{st.median(x["width"] for x in s):.1f}' if s else '-')
    print(f'{a:15} ' + ' '.join(f'{o:>24}' for o in out))
