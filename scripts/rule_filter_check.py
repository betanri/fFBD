#!/usr/bin/env python3
"""Daughter-stem rule check (validation test 6e): runs with the rule against runs without it, filtered to the samples that meet it.

A crown, stem or daughter-stem rule only removes trees; it adds no density. So the posterior with the rule must equal the posterior
without it, conditioned on the trees that meet the rule. The script runs ffbd on both inputs with several seeds, reads every sampled
tree, keeps the unconstrained samples whose fossils meet the rules, and compares the mean of every logged quantity between the two
sets (difference in standard errors; batch means, 20 batches per chain).

  rule_filter_check.py FFBD RULE_INPUT FREE_INPUT SEEDS GENS [ffbd options]
RULE_INPUT must have the daughter_stem section; FREE_INPUT is the same analysis without it. Standard library only.
"""
import sys, os, csv, subprocess, math, statistics as st

def read_rules(path):
    L = open(path).read().split("\n"); i = 0; names = {}; ext = {}; rules = []; node = {}
    while i < len(L):
        t = L[i].split()
        if t and t[0] == "tips":
            for j in range(int(t[1])): a = L[i + 1 + j].split(); names[int(a[0])] = a[1]
        if t and t[0] == "extant_nodes":
            for j in range(int(t[1])): a = L[i + 1 + j].split(); ext[int(a[0])] = a[1]
        if t and t[0] == "initial_tree":
            for j in range(int(t[1])): a = L[i + 1 + j].split(); node[int(a[0])] = (int(a[2]), int(a[3]))
        if t and t[0] == "daughter_stem":
            for j in range(int(t[1])): a = L[i + 1 + j].split(); rules.append((int(a[0]), int(a[1])))
        i += 1
    n = int([l.split()[1] for l in L if l.startswith("n ")][0])
    def ext_tips(v):
        if v < len(names): return {names[v]} if v < n else set()
        return ext_tips(node[v][0]) | ext_tips(node[v][1])
    out = []
    for f, v in rules:
        X = frozenset(ext_tips(v)); D = []
        for c in node[v]:
            u = c
            while u >= len(names) and not (ext_tips(node[u][0]) and ext_tips(node[u][1])):
                u = node[u][0] if ext_tips(node[u][0]) else node[u][1]
            D.append(frozenset(ext_tips(u)))
        out.append((names[f], X, D[0], D[1]))
    return out, set(names[i] for i in range(n))

def parse(nw):
    """Newick -> (parent dict, children dict, leaf names); internal nodes numbered."""
    par, kids, stack, cur, i, k, name = {}, {}, [], None, 0, 0, ""
    root = ("i", 0); kids[root] = []; cur = root; k = 1
    s = nw.strip().rstrip(";")
    assert s[0] == "("; i = 1; stack = [root]
    while i < len(s):
        ch = s[i]
        if ch == "(":
            nd = ("i", k); k += 1; kids[nd] = []; par[nd] = stack[-1]; kids[stack[-1]].append(nd); stack.append(nd); i += 1
        elif ch in ",)":
            if ch == ")": stack.pop()
            i += 1
        elif ch == ":":
            i += 1
            while i < len(s) and s[i] not in ",)": i += 1
        else:
            j = i
            while s[j] not in ":,)": j += 1
            nd = ("l", s[i:j]); par[nd] = stack[-1]; kids[stack[-1]].append(nd); i = j
    return par, kids

def meets(par, kids, ext, rule):
    f, X, A, B = rule
    memo = {}
    def E(u):
        if u in memo: return memo[u]
        r = frozenset([u[1]]) & ext if u[0] == "l" else frozenset().union(*[E(c) for c in kids[u]])
        memo[u] = r; return r
    def is_node(u, S):  # u is the extant node of the set S: its extant tips are S and both children hold some (or a single species)
        if E(u) != S: return False
        if len(S) == 1: return u[0] == "l"
        return all(E(c) for c in kids[u])
    anc = []; u = ("l", f)
    while u in par: u = par[u]; anc.append(u)
    return any(is_node(a, X) for a in anc) and not any(is_node(a, A) for a in anc) and not any(is_node(a, B) for a in anc)

def run(ffbd, inp, pre, seed, gens, extra):
    subprocess.run([ffbd, "run", inp, "--prefix", pre, "--seed", str(seed), "--gens", str(gens), "--print", str(gens // 10000),
                    "--tree-print", str(gens // 10000)] + extra, check=True, capture_output=True)
    lg = list(csv.DictReader(open(pre + ".log"), delimiter="\t"))
    tr = [l.split("\t")[-1] for l in open(pre + ".trees").read().splitlines()[1:]]
    return lg, tr

if __name__ == "__main__":
    ffbd, rin, fin, seeds, gens = sys.argv[1:6]; extra = sys.argv[6:]
    rules, ext = read_rules(rin); seeds = [int(x) for x in seeds.split(",")]; gens = int(gens)
    od = "rule_check_out"; os.makedirs(od, exist_ok=True)
    R, F = [], []; met_r = met_f = tot_r = tot_f = 0
    for s in seeds:
        lg, tr = run(ffbd, rin, f"{od}/rule_{s}", s, gens, extra); ch = []
        for row, t in zip(lg[1:], tr[1:]):
            p, k = parse(t); ok = all(meets(p, k, ext, r) for r in rules); tot_r += 1; met_r += ok; ch.append(row)
        R.append(ch)
        lg, tr = run(ffbd, fin, f"{od}/free_{s}", s, gens, extra); ch = []
        for row, t in zip(lg[1:], tr[1:]):
            p, k = parse(t); tot_f += 1
            if all(meets(p, k, ext, r) for r in rules): met_f += 1; ch.append(row)
        F.append(ch)
    print(f"{len(seeds)} seeds each; rules met in {met_r} of {tot_r} samples with the rule; {met_f} of {tot_f} unconstrained samples kept")
    cols = [c for c in R[0][0] if c not in ("Iteration", "Posterior", "Likelihood", "Prior", "fbd_tree")]
    def mean_se(chains, c):
        bm = []
        for ch in chains:
            x = [float(r[c]) for r in ch]; b = max(1, len(x) // 20)
            bm += [st.mean(x[i:i + b]) for i in range(0, len(x) - b + 1, b)]
        return st.mean(bm), st.stdev(bm) / math.sqrt(len(bm))
    print(f"{'quantity':24s} {'with rule':>11s} {'free, kept':>11s} {'diff in SEs':>12s}")
    worst = 0
    for c in cols:
        a, sa = mean_se(R, c); b, sb = mean_se(F, c)
        z = (a - b) / math.sqrt(sa * sa + sb * sb) if sa + sb > 0 else 0.0; worst = max(worst, abs(z))
        print(f"{c:24s} {a:11.4f} {b:11.4f} {z:12.2f}")
    print(f"largest |difference| {worst:.2f} SE over {len(cols)} quantities")
