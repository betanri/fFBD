#!/usr/bin/env python3
"""Codon-position rate variation for test 5f of VALIDATION.md (partitioned branch lengths and the young-node bias).

Real protein-coding genes differ more between codon positions (third positions are fast and saturate) than between genes.
The genes of sim_genes.py have no codon structure, so there "partition by gene" was the true scheme. Here the 23 "few" genes of
REPDIR/genes.tsv (same kinds, gene rate multipliers and base models) are re-simulated codon by codon:
  - each gene is n codons (its length rounded down to a multiple of 3); positions 1, 2 and 3 are simulated separately on the
    true phylogram scaled by gene multiplier x position multiplier, then interleaved into codons
  - position multipliers (mean 1 over the three positions), with gene-level jitter lognormal(0, 0.25) on each:
        nuclear  0.55 / 0.30 / 2.15      mitochondrial-like  0.35 / 0.12 / 2.53
  - position 3 has its own composition (GC-rich nuclear, AT-rich mito), a higher transition bias and its own gamma shape;
    position 2 has a lower gamma shape than position 1
So rates differ by gene and, more strongly, by codon position: the true scheme is gene x codon position.

  sim_codons.py REPDIR SEED OUTDIR   writes OUTDIR/codon_parts.tsv and OUTDIR/part_trees/<gene>_<pos>.nwk
Then each part is simulated with AliSim (IQ-TREE 2.3.6), and each part with AliSim, interleaves the codons and writes the partition files.
"""
import os, random, re, sys

D, seed, O = sys.argv[1], int(sys.argv[2]), sys.argv[3]
rng = random.Random(seed)
os.makedirs(f'{O}/part_trees', exist_ok=True)
tree = open(f'{D}/phylo_subs.nwk').read().strip()
rows = [l.split('\t') for l in open(f'{D}/genes.tsv').read().splitlines()[1:]]
rows = [r for r in rows if r[1] == 'few']
assert len(rows) == 23, len(rows)

POS = {'nuclear': (0.55, 0.30, 2.15), 'mito': (0.35, 0.12, 2.53)}


def scaled(t, m):
    return re.sub(r':([0-9.eE+-]+)', lambda x: f':{float(x.group(1)) * m:.10f}', t)


def freqs(base):
    f = [rng.gammavariate(100 * b, 1) for b in base]
    f = [round(x / sum(f), 4) for x in f]
    f[3] = round(1 - sum(f[:3]), 4)
    return f


with open(f'{O}/codon_parts.tsv', 'w') as fh:
    fh.write('gene\tkind\tcodons\tpos\trate_multiplier\tmodel\n')
    for g, _, kind, length, mult, alpha, model in rows:
        n = int(length) // 3
        mult, alpha = float(mult), float(alpha)
        gtr = re.search(r'GTR\{([^}]*)\}', model).group(1).split(',')
        gtr = [float(x) for x in gtr]
        pm = [p * rng.lognormvariate(0, 0.25) for p in POS[kind]]
        s = sum(pm) / 3
        pm = [p / s for p in pm]
        for k in (1, 2, 3):
            if k < 3:
                ex = [round(x * rng.lognormvariate(0, 0.2), 4) for x in gtr]
                fr = re.search(r'\+F\{([^}]*)\}', model).group(1)
                a = alpha if k == 1 else round(alpha * 0.7, 3)
            else:
                ts = 2.5 if kind == 'nuclear' else 2.0
                ex = [round(x * (ts if i in (1, 4) else 1.0) * rng.lognormvariate(0, 0.2), 4) for i, x in enumerate(gtr)]
                fr = ','.join(map(str, freqs((0.2, 0.3, 0.3, 0.2) if kind == 'nuclear' else (0.38, 0.30, 0.04, 0.28))))
                a = round(min(2.0, alpha * 2.0), 3) if kind == 'nuclear' else 0.5
            m = f'GTR{{{",".join(map(str, ex[:5]))}}}+F{{{fr}}}+G4{{{a}}}'
            r = mult * pm[k - 1]
            fh.write(f'{g}\t{kind}\t{n}\t{k}\t{r:.5f}\t{m}\n')
            open(f'{O}/part_trees/{g}_{k}.nwk', 'w').write(scaled(tree, r) + '\n')
print(f'{O}: 23 genes x 3 codon positions')
