# Fast FBD (fFBD)

A fast Bayesian program for dating phylogenies with fossils, using the **fossilized birth-death (FBD)** process.

> **Status: research software.** The math and sampling have passed strict tests, and fFBD (given a phylogram) matches full-DNA FBD (given the DNA) on two test clades. See [Validation](#validation) before using results.

## What it does

FBD dating treats fossils as real tips on the tree, instead of turning them into minimum ages on nodes. That uses all the fossils, not just the oldest ones, but standard programs become extremely slow on large trees.

fFBD (program `ffbd`) makes this practical with two ideas, and adds a third step. Two things that full-DNA FBD rarely manages on large trees follow: using all the genes, and choosing the clock by fit.

1. **The phylogram shortcut.** The DNA alignment is summarized once, as a tree of branch lengths with a standard error for each branch. The program then asks how well each proposed set of dates and rates explains those branch lengths, instead of recalculating the full DNA likelihood. Because the DNA is read only once, the branch lengths can come from all the genes of a genome-scale data set. Full-DNA FBD is in practice limited to a few dozen genes.
2. **Local updates.** When one age changes, only the affected terms are recalculated. Cost per step does not grow with tree size.
3. **Clock fitting.** Five clock models are built in. Their marginal likelihoods can be estimated on the same phylogram (`--ss`), so the clock is chosen by fit on the tree being dated, not by habit. With full-DNA FBD, this kind of fitting is practical only for small clades.

On one CPU core it runs about **2.6 million MCMC proposals per second**, with no slowdown between 1,185 and 2,424 tips.

## What is fixed and what is estimated

| Fixed | Estimated |
|---|---|
| The relationships among living species | Ages of all nodes |
| Sampling fraction of living species (rho) | Speciation, extinction and fossil sampling rates |
| Fossil clade constraints and age windows | Where each fossil attaches, within its constraint |
| Clock model (chosen by fit beforehand; see [Choosing the clock](#choosing-the-clock)) | Clock parameters (mean rate, and the rate variation of the chosen clock) |
| | Each fossil's age, within its window |
| | Whether a fossil is a sampled ancestor |
| | Rates of molecular evolution, per branch |

**Use it when** you have many fossils and a tree of living species estimated from genomic data, and full-DNA FBD is too slow.

**Tree uncertainty.** fFBD does not estimate the topology of living species: it dates a given tree. To carry topological uncertainty into the dates, date several trees and report the spread of ages across them. Good choices are trees estimated from different subsets of the genomic data (for example, different sets of loci), each by both concatenation and coalescent (species-tree) methods (as in Rincon-Sandoval et al. 2020, 2024; Santaquiteria et al. 2021, 2025). In one way this is a stronger test than dating trees from one Bayesian posterior: posterior trees all come from the same alignment and model, while trees from different subsets rest on different data. Coalescent trees give branch lengths in coalescent units, so first estimate substitution branch lengths on each topology from the alignment, for example ML branch lengths on the fixed tree; fFBD needs substitutions per site. Each tree takes hours with fFBD, so dating tens of trees is practical.

**Be careful with** very short or poorly supported branches, where a Normal approximation of branch length is doubtful.

## Install

You need a C++17 compiler and nothing else.

```bash
git clone https://github.com/betanri/fFBD.git
cd fFBD
make
```

To check the build, run `make test` (about a minute). It runs every clock on both examples with a full recomputation of the posterior after every proposal, a short stepping-stone run, and the bears example, and reports `all tests passed`.

## Quick start

The quick start uses real data: the 8 living and 14 fossil bears of the RevBayes FBD tutorial (fossils after Heath et al. 2014; branch lengths and standard errors from the tutorial's cytochrome b alignment), in `examples/bears_input.txt`. Every step takes seconds. A tiny simulated set (4 living species, 3 fossils) used by the tests is in `examples/example_input.txt`.

**Step 1. Choose the clock (highly recommended).** Fit the five clock models and list them from best to worst fit ([Choosing the clock](#choosing-the-clock) explains how to read the result):

```bash
mkdir -p runs
for c in exp ln gam igr acln; do
  ./ffbd run examples/bears_input.txt --prefix runs/clock_$c --clock $c --seed 1 --burnin 100000 \
      --ss 32 --ss-burn 100000 --ss-gens 1000000 --ss-every 200
done
grep -hv '^clock' runs/clock_*.ml.tsv | sort -t$'\t' -k7,7gr | cut -f1,7,8   # clock, stepping-stone, path sampling; best first
```

On the bears, `igr`, `acln` and `ln` tie (within 0.5 log units over three seeds) and beat `gam` by about 2.5 and `exp` by about 3.5 log units. Those gaps are below 5, so the cytochrome b data only weakly prefer the top clocks: date with the top one and check the runner-up.

**Step 2. Date with the best clock.** Run two chains:

```bash
for s in 1 2; do
  ./ffbd run examples/bears_input.txt --prefix runs/bears_r$s --clock igr --seed $s --gens 100000000 --print 50000 --tree-print 500000 --burnin 1000000
done
```

Each chain takes about 15 seconds and writes:

| File | Contents |
|---|---|
| `runs/bears_r1.log` | One row per sample: posterior, rates, and the age of each node (`age_ursinae`, ...), readable in Tracer |
| `runs/bears_r1.trees` | Sampled trees with fossils as tips |

Always run at least two chains with different `--seed` values and check that they agree. Here the two bear chains agree on node ages to about 1% (the crown of living bears at about 21 Ma).

### Choosing the clock

**Highly recommended before dating: choose the clock by fit.** fFBD has five clock models (`--clock`), and on the same DNA a different clock moved some dates by about 10 million years (test 5b), more than the phylogram shortcut ever does. Fit all five on your own phylogram, then date with the best one.

| Clock | Each branch rate |
|---|---|
| `exp` | Exponential around the mean rate |
| `ln` | Lognormal around the mean rate (uncorrelated) |
| `gam` | Gamma around the mean rate |
| `igr` | Gamma whose spread shrinks on longer branches (independent gamma rates, as in MrBayes) |
| `acln` | Lognormal around the rate of its parent branch (autocorrelated) |

To choose one, estimate the marginal likelihood of each by stepping-stone sampling (`--ss K`, with K+1 powers between posterior and prior), then keep the highest. Step 1 of the quick start shows the commands; on the example this takes a few seconds per clock.

- **Reading the result.** A difference of more than about 5 log units is strong evidence for the better clock. Smaller differences mean the data hardly prefer one clock; report the dates under the top two clocks.
- **Settings for real trees.** Use many more proposals per power: test 5c uses 64 powers, each with 250 million proposals after 50 million of burn-in. Repeat with a second seed to check that the estimate is stable.

## What you need

| Input | Notes |
|---|---|
| Rooted tree of living species | Fixed topology |
| Branch lengths and their standard errors | ML branch lengths on the fixed tree from a **partitioned** model, with bootstrap standard errors (see below) |
| Fossil list | Name, minimum and maximum age |
| Clade constraints | Which fossils belong in which clades |
| Crown and stem rules (optional) | Where a fossil may sit relative to a crown group: four rules (see below) |
| Node calibrations (optional) | Hard minimum and maximum ages for chosen nodes, the root included (see below) |
| Sampling fraction (rho) | Share of living species included |
| Root age prior | `uniform MIN MAX`, `soft MIN MAX SD P` (soft bounds) or `normal MEAN SD` |
| Clock model | Chosen by fitting the five clocks on your phylogram (highly recommended; [Choosing the clock](#choosing-the-clock)), then set with `--clock` |
| Starting tree (optional) | **Built automatically** by `scripts/ffbd_convert.py`. Give your own `initial.nwk` only if you want a particular start |

**Partition the branch lengths.** Fit them by gene and codon position, or with a scheme merged by ModelFinder (IQ-TREE `-p partitions.nex -m MFP+MERGE`). A single model for all sites gets short branches wrong: young nodes come out too old relative to old ones, and their intervals miss the true age more often (test 5f).

These go in one text file. `examples/example_input.txt` shows every section: the settings, then `tips`, `initial_tree`, `constraints`, `extant_nodes`, `branches` and `terms`. If your inputs are already in RevBayes-style files, `scripts/ffbd_convert.py` builds this file for you.

**You do not need a starting tree.** MCMC needs one valid tree to start from: fossils attached, every constraint, crown or stem rule and calibration satisfied. If the input folder has no `initial.nwk`, the converter builds one:

- each fossil goes on the stem of its clade, at its minimum age;
- crown-rule and daughter-stem fossils go inside the crown, on the stem of its first daughter group;
- node ages are pushed up wherever a fossil or calibration needs it.

It prints `starting tree built automatically`. `--write-start FILE` saves the tree, and `--auto-start` uses it even when `initial.nwk` exists. The start does not affect the result: from an automatic and from a supplied start, the node ages agree within chain noise (test 6d).

### Crown and stem rules

Every fossil follows one of four rules. The clade constraint alone gives the first; three optional sections narrow it, each for a node X of the living tree:

| Rule | Where the fossil may attach | Input section | Converter file |
|---|---|---|---|
| **Total group** (default) | Anywhere in its clade: on the stem or inside the crown | none (the clade constraint) | none |
| **Stem** | On the stem of X only, never inside crown X | `stem` | `stem.tsv` |
| **Crown** | Anywhere inside crown X | `crown` | `crown.tsv` |
| **Daughter stem** | Inside crown X but in neither daughter crown: on the stem of one of X's two daughter groups | `daughter_stem` | `daughter_stem.tsv` |

- **Stem.** The clade constraint keeps the fossil in the total group, so it sits on the stem.
- **Crown.** Crown X must then be at least as old as the fossil.
- **Daughter stem.** Use it when a fossil has the features of crown X but none of either daughter group, so it sits near the first split inside X without saying on which side. The daughter groups are the first nodes of the living tree below X on each side (a single living species counts). A fossil may also sit on a side branch of another fossil that hangs from those stems. The rule holds the fossil close to node X, so it acts much like a node calibration: use it only where the evidence rules out both daughter crowns.

Each section is the section name and K, followed by K lines `fossil node`. `fossil` is the tip index and `node` the extant node index (X) used in the input file. The converter reads the files above (fossil name, extant node name) from the input folder. A fossil can have only one rule. The rules only restrict where a fossil can attach; they add no age prior. `examples/example_crown_stem.txt` shows the crown and stem rules, and `examples/bears_daughter_stem.txt` the daughter-stem rule.

### Node calibrations

A node of the living tree, the root included, can be held between a minimum and a maximum age. Add a section `calibrations K` followed by K lines `node min max` (extant node index, ages in Ma). The bounds are hard: the node never leaves them, and inside them its age follows the FBD density and the data, with a uniform calibration density. The converter reads them from `calibrations.Rev` in the input folder, one line per node in the RevBayes form `cal_X ~ dnUniform(age_X - MAX, age_X - MIN)` (X = extant node name).

Fossils already date the tree, so most analyses need no calibrations. One case where they help is the root: if the root has few fossils close to it, its age can swing widely between chains or runs. A hard maximum on the root (with the root prior kept) stops that.

## Main options

| Option | Meaning |
|---|---|
| `--prefix P` | Start of the output file names |
| `--seed S` | Random seed |
| `--gens G` | Number of proposals |
| `--burnin N` | Proposals used to tune moves, not logged |
| `--print N`, `--tree-print N` | How often to log parameters and trees |
| `--checkpoint N`, `--resume` | Save progress and continue a run exactly where it stopped |
| `--mode control` | Run without molecular data, to see what the fossils and priors alone imply |
| `--log-rates` | Also log branch rates |
| `--clock exp\|ln\|gam\|igr\|acln` | Clock model. Always set it, after fitting the clocks ([Choosing the clock](#choosing-the-clock)). If omitted, `exp` is used and a warning is printed |
| `--ss K --ss-gens G` | Instead of a normal run, estimate the marginal likelihood of the phylogram by stepping-stone sampling with K+1 powers, to compare clock models |
| `--ss-power fbd` | With `--ss`: put the powers on the FBD density instead, to compare rate models (no shifts, fixed or estimated shifts) |
| `--rate-shifts SPEC` | Rate shifts: fixed times, `free:K`, windows `lo:hi`, or `mass-extinctions` (The skyline model, below) |
| `--heated K` | Run K heated copies of the chain and log the cold one, for posteriors with several peaks (Heated chains, below) |
| `--check` | Check internal consistency while running (slow) |

## The model

- Fossilized birth-death process with sampled ancestors, conditioned on the root age (or the origin). Rates are constant by default, or piecewise-constant between shift times, fixed or estimated (skyline FBD, below). The density matches RevBayes 1.3.2 `dnBDSTP` (with `timeline` when there are shifts) inside `dnConstrainedTopology`.
- Fossil ages are uniform within their windows.
- Data: for each branch, observed length ~ Normal(rate x time, standard error). The two branches at the root share one observation.

### The clock model (rate changes across branches)

The molecular substitution rate can differ from branch to branch, so lineages and clades can evolve at different speeds. There is no built-in trend through time.

- One rate per branch of the living-species tree, with mean rate m ~ Exp(10). Five clock priors (`--clock`):
    - `exp`: each rate ~ Exponential(mean m).
    - `ln`: uncorrelated lognormal, mean m.
    - `gam`: each rate ~ Gamma(shape k, mean m).
    - `igr`: independent gamma rates as in MrBayes; the rate variance shrinks on longer branches.
    - `acln`: autocorrelated lognormal; each branch's rate is centered on its parent branch's rate, with variance growing with branch duration.
- The clock with the best fit can be chosen by marginal likelihood (`--ss`; test 5c).

### The skyline model (rate changes through time)

Speciation, extinction and fossil sampling can change through time at shift times, fixed or estimated from the data. At any time they are the same for all lineages alive, so they do not differ across clades. This is skyline FBD, as in BEAST2, RevBayes `dnBDSTP(timeline=...)` and MrBayes 3.2. Add one line to the input (or use `--rate-shifts` on the command line, which overrides it):

```
rate_shifts 66 201.3
```

- K times in Ma, strictly increasing, give K+1 intervals; interval 1 is the youngest. Each interval has its own lambda, mu and psi.
- `shift_rates lambda,mu` lets only those rates change; a rate left out (here psi) is shared by all intervals.
- Each free rate ~ Exponential(rate 10), as in the constant model, with no smoothing between intervals: an interval with few nodes or fossils stays close to the prior.
- rho applies at the present only; there are no mass extinctions or extra sampling at the shifts.
- The log gets one column per interval: `speciation_rate[1]`, `extinction_rate[1]`, `psi[1]`, ... (youngest first).
- Without shifts the output is byte-identical to the constant-rate model. With 2 shifts a run takes about 60% longer per proposal.
- Validation: test 1b of VALIDATION.md (density equals RevBayes on 1,012 states, cache check, prior recovery, simulation-based calibration).

**Estimated shift times.** The shift times can be parameters of the model:

| `--rate-shifts` | Shift times |
|---|---|
| `66,201.3` | Fixed (as above) |
| `free:K` | K times estimated anywhere in a window, ordered; prior uniform (density K!/W^K). Window: `--shift-window lo,hi`, default 0 to the oldest fossil's minimum age |
| `60:72,195:207` | One time estimated in each window, prior uniform on the window (windows must not overlap). A plain number in the list stays fixed |
| `mass-extinctions` | Fixed at the major extinction events younger than the oldest fossil: K-Pg 66.0, end-Triassic 201.4, end-Permian 251.902, end-Devonian 358.9, Late Devonian (Frasnian-Famennian) 371.1, end-Ordovician 443.8 Ma (ICS chart) |

- Three moves change the times. `shift_time` slides one time or redraws it between its neighbours. `shift_data` proposes a time where the fossils and the tree make a change in fossil sampling likely, and can redraw the fossil-sampling rates of the two intervals with it. `shift_relocate` (`free:K`, K of 2 or more) takes one shift out and puts it back anywhere in the window, so a shift can jump over the others. Each move has weight 5 per estimated time. All are exact Metropolis-Hastings moves. The log gets `shift_time[1]`, ... (youngest first).
- When fossil sampling may change at the shifts, an estimated shift tends to land on a rich fossil bed, where many fossils date to a short interval, rather than on a change in speciation or extinction. With `shift_rates lambda,mu` fossil sampling stays the same in every interval, so the shifts are about speciation and extinction only.
- Each extra shift adds an interval with three more rates. Choose the number of shifts by marginal likelihood (below), not by adding shifts until the dates stop changing.
- Validation: tests 1c and 2d of VALIDATION.md (exact numerical integration, simulation-based calibration, recovery of the number of shifts).

**Choosing the rate model.** Rate models (no shifts, K estimated shifts, fixed or preset shifts) are compared by the marginal likelihood of the whole model, including the FBD density of the tree and fossils. Use `--ss K --ss-power fbd`: the stepping-stone powers go on the FBD density, while the phylogram stays at full weight. Everything else has the same prior in every rate model, so differences of `stepping_stone` in `PREFIX.ml.tsv` are log Bayes factors between rate models. (The default `--ss-power data` puts the powers on the phylogram; that compares clocks, not rate models.) A difference above 3 log units is strong support, above 5 very strong (Kass and Raftery 1995). Run two replicates of each model; they should agree within about 1 log unit, or the run needs more powers or more proposals per power.

```
for m in none free:1 free:2 free:3 mass-extinctions; do
  ffbd run input.txt --prefix ss_${m/:/} --clock gam --burnin 1000000 --rate-shifts $m --ss 100 --ss-power fbd --ss-gens 2000000 --ss-burn 200000 --ss-every 1000
done
```

On large trees the FBD density changes steeply near full power, so they need many more powers and longer runs per power than small ones (several hundred powers on thousands of tips). Check that the replicates agree before comparing models.

### Heated chains

On a large tree the posterior of estimated shift times can have several separate peaks, and a single chain can stay on one of them. Two chains with different seeds that disagree on the shift times show this. `--heated K` runs K copies of the chain side by side (Metropolis-coupled MCMC). Copy r targets the posterior with the FBD density raised to the power H^(r/(K-1)), H = `--hottest` (default 0.3), so the hotter copies cross between peaks more easily. Every `--swap-every` proposals two neighbouring copies may swap states. Only the cold copy (power 1) is logged, so the output is a sample of the posterior.

| Option | Meaning |
|---|---|
| `--heated K` | Number of copies (2 or more) |
| `--hottest H` | Power of the hottest copy, between 0 and 1 (default 0.3) |
| `--heat-on fbd\|all` | Heat the FBD density only (default), or the FBD density and the phylogram |
| `--swap-every N` | Proposals of each copy between swap attempts (default 1000) |
| `--threads T` | Threads (default one per copy); the output does not depend on it |

- The end of the run prints the swap acceptance between neighbouring copies. If one pair rarely swaps (below about 0.05), the heats are too far apart: add copies or raise `--hottest`. On a tree of thousands of tips the FBD density is very sharp, so neighbouring copies need powers that differ by about 0.01.
- Heated chains also work with `--ss` (the stepping-stone power multiplies the heats) and with `--checkpoint` and `--resume`.
- A run with K copies does K times the work. Use it when two plain chains disagree, and still run two heated chains with different seeds and check that they agree.
- Validation: test 2e of VALIDATION.md.

## Validation

All numbers, per-test details and result files are in [VALIDATION.md](VALIDATION.md) and `validation/`, with the same test numbers as here.

fFBD is new code, and it replaces the DNA with a phylogram. Both need to earn trust. The validation follows the standard recipe for Bayesian software (Cook et al. 2006; Talts et al. 2018; Mendes et al. 2025) and for FBD simulation studies (Heath et al. 2014; Barido-Sottani et al. 2019). It is organized as six questions. For each one: why it matters, how it was tested, and the answer.

### Short answers

| # | Question | Answer |
|---|---|---|
| 1 | Does fFBD compute the FBD probability correctly? | **Yes.** Every term matches RevBayes to within 4e-9 (rounding error). |
| 2 | When the true answer is known, does it find it? | **Yes.** Simulations with known truth pass on trees of 2 to 8,056 living species. |
| 3 | Does an independent program give the same answer? | **Yes.** MrBayes agrees within 0.5%. RevBayes 1.3.2 differs because of one move; with that move off, it agrees too. |
| 4 | Does the phylogram shortcut give the same dates as the DNA? | **Yes, under the same model.** The differences are as small as those between DNA chains of the same analysis. |
| 5 | What changes the dates more than the shortcut does? | **The genes behind the branch lengths, and the clock model.** The sampling fraction matters less. |
| 6 | Does it work on large trees? | **Yes.** Chains converge at 2,424 species, and speed does not drop with tree size. |

### Three kinds of data

fFBD never reads DNA. Its data are a **phylogram**: the branch lengths of the fixed tree of living species, in substitutions per site, each with a standard error, estimated once from the DNA. Each test is labeled by what it uses:

- **No data**: fossils and priors only. This isolates the tree and fossil part of the model, where programs differ most.
- **Phylogram**: branch lengths and standard errors, estimated from DNA or simulated.
- **DNA**: the full alignment, analyzed by RevBayes or MrBayes as the reference. fFBD then gets a phylogram estimated from that same DNA.

The real test clades are Gobiaria (24 living species, 10 fossils) and Protacanthopterygii (22 living species, 14 fossils), each with 23 genes (16,107 sites), plus catfishes (713 living species, 94 fossils) and the RevBayes bear tutorial.

### 1. Does fFBD compute the FBD probability correctly?

**Why it matters.** Every result rests on the probability of a dated tree given the fossils and branch lengths. A small error there biases everything.

**How it was tested.** *No data and phylogram.* Fixed trees and parameter values were scored by fFBD and by RevBayes, the established FBD program, term by term: tree density, root prior, rate priors, branch-length terms and fossil windows. 89 states on 4 data sets.

**Answer.** The totals match to within 4e-9, which is rounding error (test 1). With rate shifts the density matches RevBayes too (test 1b). With estimated shift times, the posterior of the times and the marginal likelihood match exact numerical integration (test 1c).

### 2. When the true answer is known, does it find it?

**Why it matters.** Correct probabilities are not enough; the sampler also has to explore them correctly. Bugs in proposals or bookkeeping show up only here.

**How it was tested.** *No data, and simulated phylograms.* Hundreds of data sets were simulated from the model with known true values, using separate simulation code. If the sampler is correct, the true value falls at a random position among the posterior samples, so its rank is uniform. This is simulation-based calibration (SBC).

- **2a.** Small trees, 500 replicates in each of 3 settings, 13 quantities.
- **2b.** Realistic trees, 500 replicates of 2 to 8,056 living species and up to 15,838 fossils, 9 quantities.
- **2c.** A 3-species case with an exact answer.
- **2d.** Rate shifts: simulation-based calibration with estimated shift times (6 settings, 1,000 replicates each), and recovery of the number of shifts by marginal likelihood on 90 simulated trees. Repeated with the new shift-time moves.
- **2e.** Heated chains: simulation-based calibration (4 settings, 1,000 replicates each, then 2,000 more in each of the 2 settings with the most failures), and heated against plain chains on one data set (20 runs each).

**Answer.** All passed. In 4 of the first 9 heated-chain runs one quantity fell below the threshold; a larger repeat with new seeds and the direct comparison with plain chains found no difference (VALIDATION.md, test 2e). Testing found and fixed two bugs along the way (a cached rate sum and a checkpoint restart), and all tests were rerun after the fixes.

### 3. Does an independent program give the same answer?

**Why it matters.** Two programs written from the same equations can still share a mistake in the way they move through tree space. Agreement between independent code is strong evidence.

**How it was tested.** *No data.* fFBD, MrBayes 3.2.7 and RevBayes 1.3.2 got the same fixed tree of living species, the same fossils with age windows, and the same fixed speciation, extinction, fossil sampling and root age (test 3a). With no data, all three should give the same node ages and the same number of sampled ancestors, meaning fossils that sit directly on a branch rather than on a side branch. Test 3b reran RevBayes with one move switched off.

**Answer.** fFBD and MrBayes agree (node ages within 0.4 to 0.6%, the same number of sampled ancestors). RevBayes finds about one extra sampled ancestor and makes young nodes up to a third older. The likely cause is the RevBayes move `mvCollapseExpandFossilBranch`, which is not reversible when fossils have age ranges (details in VALIDATION.md). With it off, RevBayes agrees with fFBD within 0.4% (test 3b). The RevBayes developers have not yet confirmed this. BEAST2 was not used, because it cannot hold the tree of living species fixed while fossils move.

![fFBD, MrBayes and RevBayes on the same prior](validation/three_programs_no_data.png)

*Test 3a.* Each dot is a node. Top: how far each program is from MrBayes, in %. Blue (fFBD) sits on zero; orange (RevBayes) sits above it, mostly on young nodes. Bottom: the number of sampled ancestors per tree. fFBD and MrBayes match; RevBayes is shifted by about one.

### 4. Does the phylogram shortcut give the same dates as the DNA?

**Why it matters.** This is the question that matters most in practice. The shortcut is what makes fFBD fast, so it must not change the answer.

**How it was tested.**

- **4a. Same DNA, same model.** *DNA vs phylogram.* RevBayes ran full-DNA FBD on the 23 genes of each test clade (6 chains of 22 to 25 hours each, with the move from test 3 removed). fFBD got the phylogram estimated from the same 23 genes, with the same fossils, priors and clock.
- **4b. Simulated DNA, known true ages.** *Phylogram from simulated DNA.* 100 simulated trees of 30 to 120 species, with DNA simulated on them, analyzed exactly as for real clades.
- **4c. Catfishes.** *DNA vs phylogram.* An earlier full-DNA analysis of 713 species. It is not like for like, because it used another clock and the RevBayes move from test 3.

**Answer.**

- **4a.** The ages agree to 1.3 and 1.8% at the median, which is no more than the DNA chains differ from each other (1.5 and 3.2%). Roots agree to 0.1 Myr.
- **4b.** The true ages fall inside the 95% interval 91% of the time with the shortcut, against 93% with exact branch lengths. Nodes under 5 Ma come out about 1 Myr too old with or without the shortcut, so that bias comes from the model and priors.
- **4c.** The ages are close (r = 0.994), but not identical.

![Node ages, full-DNA FBD against fFBD](validation/node_ages_full_vs_fast.png)

*Tests 4a and 4c.* Each dot is a node. Horizontal: its age from full-DNA FBD. Vertical: its age from fFBD. On the diagonal means the same age.

### 5. What changes the dates more than the shortcut does?

**Why it matters.** If the shortcut is not the weak link, users should know what is.

**How it was tested.**

- **5a. Which genes.** *DNA vs several phylograms.* Everything stayed the same except the phylogram given to fFBD: one from the same 23 genes, one from all genes (about 1,100 genes, 580,000 sites), and one from a genome-scale phylogram of a larger tree, pruned to the clade.
- **5b. Which clock, by a second DNA program.** *DNA vs phylogram.* MrBayes ran full-DNA FBD on the same 23 genes. MrBayes cannot run the clock that fFBD and RevBayes use (one exponential rate per branch), so it used independent gamma rates (IGR).
- **5c. Which clock fits best.** Clocks are compared by marginal likelihood, which measures how well each model predicts the data while penalizing extra flexibility (stepping-stone sampling; Xie et al. 2011). There are two complementary fits:
    - *DNA.* RevBayes on the 23 genes of the two test clades, with six clocks: strict, exponential (the clock used in tests 1 to 5b), lognormal, gamma, IGR and autocorrelated lognormal. Node ages are either fixed at the DNA estimates or estimated with fossil minimum ages.
    - *Phylogram.* fFBD itself (`--ss`) on the 1,184-species phylogram used for dating, with all five fFBD clocks.

    The best-fitting clock will be used for the final dates. With full-DNA FBD, this kind of clock fitting is practical only on small clades (see "Why fFBD may do better" below); fFBD can do it on the full tree being dated. Three of the clocks (gamma, IGR, autocorrelated lognormal) and the stepping-stone code are new, so they first passed the same kinds of checks as questions 1 and 2:
    - their prior densities match an independent implementation (scipy) to rounding error;
    - with no data they recover their priors exactly;
    - on a small example, the stepping-stone estimate matches a brute-force average of the likelihood over 500,000 prior draws within Monte Carlo error;
    - with the exponential clock, the program gives byte-identical output to the previous version.
- **5d. Sampling fraction.** *Phylogram.* fFBD fixes rho, the fraction of living species in the tree. A 1,184-species phylogram (pruned from a genome-scale phylogram of a larger tree) with 1,588 fossils was dated at rho = 3.25, 10, 15, 25 and 50%.
- **5e. More genes, known truth.** *Phylograms from simulated DNA.* Gene sets were simulated on the 100 trees of test 4b. The genes differ the way real genes do: each has its own rate, base composition and rate variation, and a few are fast, mitochondrial-like genes. There are two sets: 23 genes (about 18,000 sites), and 1,000 genes (about 500,000 sites) that include the 23. Each set was analyzed like the real clades: concatenated, with ML branch lengths on the true tree.
- **5f. Partitioned branch lengths, known truth.** *Phylograms from simulated DNA.* The same trees, with genes that also differ by codon position. Branch lengths were fitted with one model for all sites, by gene, by codon position, by gene x codon position, and with a ModelFinder-merged scheme.

**Answer.**

- **5a. The genes matter.** The 23-gene phylogram reproduces the 23-gene DNA analysis (1.0 and 1.5% at the median). Phylograms from all genes make nodes 2.5 to 3.8 Myr younger on average, still inside the DNA 95% intervals. How the standard errors were computed makes no difference.
- **5b. The clock can matter.** On Gobiaria, MrBayes agrees with both RevBayes and fFBD (about 3%). On Protacanthopterygii, the two DNA programs disagree: MrBayes puts the inner nodes about 10 Myr younger, with intervals twice as wide, while fFBD stays with RevBayes, whose clock it shares.
- **5c. Gamma or lognormal, not exponential.** The table gives each clock's log marginal likelihood minus that of the best clock in its column (0 = best; more negative = worse fit; a gap above about 5 log units is strong evidence). DNA columns use ages fixed at the DNA estimates, mean of 2 replicates.

    | Clock | DNA, Gobiaria | DNA, Protacanthopterygii | Phylogram, 1,184 species |
    |---|---|---|---|
    | Gamma | -1.4 | -1.6 | **0** (two seeds: 5,746 and 5,761) |
    | Lognormal | **0** | **0** | -5.3 |
    | Autocorrelated lognormal | -14.9 | -2.7 | -42 |
    | Exponential | -8.7 | -12.2 | -147 |
    | IGR | -27.3 | -35.5 | -464 |
    | Strict | -740 | -523 | not in fFBD |

    - **The top two.** Gamma and lognormal lead in every column, within about 5 log units of each other, which is about the gap between the two gamma seeds. The data cannot tell them apart.
    - **The rest.** Exponential, the clock used in tests 1 to 5b, is clearly worse everywhere, decisively on the large phylogram. IGR, the MrBayes clock of test 5b, fits worst of the relaxed clocks. A strict clock is far behind.
    - **Choice.** The final dates use gamma, with lognormal as a check.
    - **Run time.** Each phylogram fit took 2.5 to 3.8 hours on one core. A DNA fit of a single clock on 22 to 24 species took 1.5 to 10 hours on 8 cores in RevBayes.
    - **Cross-check.** With node ages fixed at the MrBayes estimates, or estimated under fossil minimum ages, the relaxed clocks again beat exponential and strict. There the autocorrelated clock leads for Protacanthopterygii, with lognormal within 2.5 log units. The remaining free-age runs will be added when done.
- **5d. Rho matters less.** Raising rho from 3.25 to 50% makes crown ages younger: under 1% at the root (held by its prior) and up to 8% on younger crowns. At rho = 10%, no crown moves more than 2%.
- **5e. More genes are more accurate, modestly.** With 1,000 genes, node ages are closer to the truth on 63% of nodes. The median error drops from 3.55 to 3.39 Myr (3.13 with perfect branch lengths), and 95% coverage rises from 89.5 to 90.7% (93.1% with perfect branch lengths). The gain is largest on nodes under 5 Ma. A 23-gene phylogram gives the same answer as 23-gene DNA (tests 4a and 5a), so this compares full-DNA FBD on the genes it can afford with fFBD on genome-scale data.
- **5f. Partition the branch lengths.** One model for all sites makes short branches too long, so young nodes come out too old: under 5 Ma, coverage 0.79 and +106% bias. Partitioning by gene x codon position, or the ModelFinder scheme, gives 0.90 and +67%, close to perfect branch lengths (0.93 and +60%). By codon alone or by gene alone helps much less. Older nodes hardly change.

![Phylograms against 23-gene DNA](validation/phylogram_vs_dna_controlled.png)

*Test 5a. Does fFBD give the same dates as a full DNA analysis, and what changes them?* Each dot is one node of two fish clades. The x-axis is its age in the 23-gene DNA analysis (RevBayes); the y-axis is how far another run's age is from it, in % (0 = same age, below 0 = younger). Every run has the same tree, fossils and priors; only the input changes. **Grey:** a second set of DNA chains, i.e. the noise you get just by rerunning the DNA analysis. **Blue:** fFBD given a phylogram from the same 23 genes. It sits with grey (median 1.0 and 1.5%), so replacing the DNA with a phylogram costs nothing. **Orange and green:** fFBD given phylograms from far more data (all ~1,100 genes, or a larger genome-scale tree pruned to the clade). Their nodes are younger, by 5 to 9% at the median and most on the youngest nodes, because more genes give different branch lengths, not because fFBD differs from DNA. Test 5e checks against simulated true ages which answer is closer to the truth.

### 6. Does it work on large trees?

**Why it matters.** Large trees are the reason fFBD exists.

**How it was tested.** *Phylogram.*

- **6a. Convergence.** Two independent chains each were run at 1,185 and 2,424 species, with 20 billion proposals each (about 2 hours on one core).
- **6b. Crown and stem rules.** The optional rules that keep a fossil inside a crown, on its stem, or on the stem of one of its two daughter groups were compared with unconstrained runs filtered to the same rules.
- **6c. Speed.** Speed was measured at both tree sizes.
- **6d. Automatic starting tree.** Starting trees were built automatically for two small clades and for two 1,184-species sets (1,587 and 1,588 fossils, 430 constraints, crown and stem rules, node calibrations). The same analysis was then run from the automatic start and from a supplied one.

**Answer.** The chains converge (effective sample sizes 500 to 800 for the main parameters; the two chains agree on node ages to 0.9% at the median). The rules behave exactly as intended (all quantities within 1.4 standard errors for the crown and stem rules and 1.7 for the daughter-stem rule). Every automatic starting tree passes all checks. From the automatic and the supplied start, node ages agree to 0.17 and 0.29% at the median (largest 0.7 and 1.1%), within chain noise. Speed is the same at both sizes, about 2.6 million proposals per second.

### What the validation shows, and what it does not

**It shows that:**

- fFBD computes and samples the FBD model correctly (questions 1 and 2).
- An independent program, MrBayes, gives the same answer (question 3).
- Given the same DNA and model, the phylogram shortcut gives the same dates as analyzing the DNA (question 4).

**It does not show that:**

- **The dates are true.** The true ages of real clades are unknown. More genes and better clocks reduce error, but model error remains.
- **The clock is right.** The clock and the genes behind the branch lengths can move dates more than the shortcut does (question 5), so they deserve the same care as fossil choice.
- **Very young nodes are accurate.** Nodes under about 5 Ma tend to come out about 1 Myr too old. Part of this comes from the model and priors, and full-DNA FBD shares it; the rest comes from unpartitioned branch lengths and mostly goes away with partitioning (test 5f).
- **The topology is right.** Each run dates one tree of living species. Topological uncertainty is covered by dating several trees (see Tree uncertainty, above).

### Why fFBD may do better than full-DNA FBD on real data

Full-DNA FBD is slow, and that forces three compromises in practice. fFBD removes all three.

**1. All the genes, not a few dozen.**

- **Full-DNA FBD in practice** is run on a few dozen genes or fewer. Here, RevBayes on 23 genes and about two dozen species took 22 to 25 hours per chain and still mixed poorly, so six chains had to be pooled. Running the same analysis on 36 times more sites, or on thousands of species, is out of reach.
- **fFBD** reads the DNA only once, as branch lengths, so those can come from every gene of a genome-scale data set.
- **What test 5a shows.** Given the same genes, fFBD matches full-DNA FBD, so the shortcut does not change the answer. The data do change it: all-gene branch lengths give ages a few million years younger, and they are 5 to 6 times more precise (standard error relative to branch length).
- **What test 5e shows.** With known true ages, 1,000 genes give more accurate dates than 23: closer to the truth on 63% of nodes, with the largest gain on young nodes. The simulated genes follow the analysis model, so the gain on real data may be smaller or larger.

**2. The clock chosen by fit, not by habit.**

- **Full-DNA FBD in practice.** RevBayes and MrBayes offer several clocks and stepping-stone sampling, but on DNA fitting them is practical only for small clades. Here, one RevBayes stepping-stone run for one clock on 22 to 24 species and 23 genes took 1.5 to 10 hours on 8 cores. Repeating that for every clock on 1,184 species and 1,588 fossils is out of reach, so studies pick one clock and rarely test it.
- **fFBD** runs the whole pipeline on the phylogram of the tree actually being dated, one core per clock (2.5 to 3.8 hours each on 1,184 species and 1,588 fossils): it fits five clocks by marginal likelihood, keeps the best, and dates with it.
- **What test 5b shows.** The clock matters. On the same DNA, a different clock moved some dates by about 10 million years.
- **What test 5c shows.** Exponential, which most FBD analyses would simply have used here, fits clearly worse than gamma or lognormal on every data set: by 9 to 22 log units on the DNA of the two small clades, and by about 147 on the 1,184-species phylogram. Each phylogram fit took 2.5 to 3.8 hours on one core.

**3. Rate shifts chosen by fit, not by habit.**

- **Full-DNA FBD in practice.** BEAST2, RevBayes and MrBayes can fit skyline FBD, but comparing no shifts against one, two or three shifts, or against shifts at extinction events, needs a stepping-stone run per model on the DNA: the same cost that limits clock fitting to small clades. Studies therefore keep rates constant or place the shift times by hand. RevBayes 1.3.2 also cannot estimate the shift times themselves (test 1c).
- **fFBD** fits each rate model on the phylogram of the tree being dated (no shifts, K shifts at times estimated from the data, or shifts at preset extinction events) and keeps the one the data support (`--ss-power fbd`).
- **What test 2d shows.** On simulated trees of 30 to 120 living species, the marginal likelihood picked no shift on 29 of 30 trees simulated without a shift (30 of 30 with the new moves), and one shift on 27 of 30 simulated with one, with the true shift time inside the 95% interval on 26 of 30. It is conservative: it rarely adds a shift that is not there, and on trees this small it seldom detects a second one.

### Reproducing the tests

| Test | Script | Result file |
|---|---|---|
| 1 | `density_check.py` | `density_check_A06_A12.txt`, `density_check_bears_catfish.txt` |
| 2a | `sbc.py` | `sbc_with_data_500reps.txt`, `sbc_prior_only_500reps.txt`, `sbc_origin_mode_500reps.txt` |
| 2b | `sim_study.py sbc`, `sbc_ranks.py` | `sbc_realistic_500reps.txt`, `sbc_rank_histograms.png` |
| 1b | `sky_density_check.py`, `sbc_sky.py` | `skyline_checks.txt` |
| 1c | `est_grid_check.py`, `est_rb_check.py` (fixed shifts only) | `skyline_estimated_checks.txt` |
| 2c | `analytic_check.py` | `analytic_check.txt` |
| 2d | `sbc_sky.py --shifts free:K:lo:hi`, `sim_study.py simulate --shifts` | `skyline_estimated_checks.txt`, `shift_moves_heated_checks.txt` |
| 2e | `sbc_sky.py --ffbd-extra "--heated 3 ..."` | `shift_moves_heated_checks.txt` |
| 3a | `make_fixed3.py`, `compare_fixed3.py` | `three_program_prior_A06_A12.txt`, `compare_pilots_summary.txt` |
| 3b | `compare_runs.py` | `E1_noCE_summary.txt` |
| 4a, 4c | `branch_scatter.py --ages-only` | `node_ages_full_vs_fast.tsv` |
| 4b | `sim_study.py simulate`, `sim_study_summary.py` | `sim_study_100reps.txt`, `simulation_coverage_by_age.png` |
| 5a | `ctrl_compare.py` | `phylogram_vs_dna_controlled.txt` |
| 5b | `make_mb_dna.py`, `compare_mbdna.py` | `mrbayes_dna_A06_A12.txt` |
| 5c | `clock_ss.py` (DNA, RevBayes); `ffbd --clock C --ss K` (phylogram); `clock_density_check.py` | `clock_checks.txt`, `clock_fit_results.txt` |
| 5e | `sim_genes.py`, `sim_genes_summary.py` | `sim_genes_100reps.txt` |
| 5f | `sim_codons.py`, `young_fix_summary.py`, `codon_test_summary.py` | `young_fix_50reps.txt`, `codon_test_50reps.txt` |
| 6b | `examples/example_crown_stem.txt`, `examples/bears_daughter_stem.txt`, `rule_filter_check.py` | `crown_stem_test_11seeds.txt`, `daughter_stem_test_11seeds.txt`, `regression_previous_release.txt` |
| 6d | `ffbd_convert.py --auto-start` | `auto_start_check.txt` |
| Figures | `validation_figures.py` | `validation/*.png` |

### References

- Barido-Sottani J, Aguirre-Fernández G, Hopkins MJ, Stadler T, Warnock R (2019). Ignoring stratigraphic age uncertainty leads to erroneous estimates of species divergence times under the fossilized birth-death process. *Proceedings of the Royal Society B* 286: 20190685.
- Cook SR, Gelman A, Rubin DB (2006). Validation of software for Bayesian models using posterior quantiles. *Journal of Computational and Graphical Statistics* 15: 675-692.
- Gavryushkina A, Welch D, Stadler T, Drummond AJ (2014). Bayesian inference of sampled ancestor trees for epidemiology and fossil calibration. *PLoS Computational Biology* 10: e1003919.
- Heath TA, Huelsenbeck JP, Stadler T (2014). The fossilized birth-death process for coherent calibration of divergence-time estimates. *PNAS* 111: E2957-E2966.
- Höhna S, Landis MJ, Heath TA, et al. (2016). RevBayes: Bayesian phylogenetic inference using graphical models and an interactive model-specification language. *Systematic Biology* 65: 726-736.
- Mendes FK, Bouckaert R, Carvalho LM, Drummond AJ (2025). How to validate a Bayesian evolutionary model. *Systematic Biology*, doi:10.1093/sysbio/syae064.
- Rincon-Sandoval M, Duarte-Ribeiro E, Davis AM, et al. (2020). Evolutionary determinism and convergence associated with water-column transitions in marine fishes. *PNAS* 117: 33396-33403.
- Rincon-Sandoval M, et al. (2024). Ecological diversification of sea catfishes is accompanied by genome-wide signatures of positive selection. *Nature Communications* 15: 10040.
- Ronquist F, Teslenko M, van der Mark P, et al. (2012). MrBayes 3.2: efficient Bayesian phylogenetic inference and model choice across a large model space. *Systematic Biology* 61: 539-542.
- Santaquiteria A, Siqueira AC, Duarte-Ribeiro E, et al. (2021). Phylogenomics and historical biogeography of seahorses, dragonets, goatfishes, and allies (Teleostei: Syngnatharia): assessing factors driving uncertainty in biogeographic inferences. *Systematic Biology* 70: 1145-1162.
- Santaquiteria A, Miller EC, Rosas-Puchuri U, Pedraza-Marrón CR, Troyer EM, Westneat MW, Carnevale G, Arcila D, Betancur-R R (2025). Colonization dynamics explain the decoupling of species richness and morphological disparity in syngnatharian fishes across oceans. *The American Naturalist* 205(3), doi:10.1086/733931.
- Talts S, Betancourt M, Simpson D, Vehtari A, Gelman A (2018). Validating Bayesian inference algorithms with simulation-based calibration. arXiv:1804.06788.
- Xie W, Lewis PO, Fan Y, Kuo L, Chen MH (2011). Improving marginal likelihood estimation for Bayesian phylogenetic model selection. *Systematic Biology* 60: 150-160.
- Zhang C, Stadler T, Klopfstein S, Heath TA, Ronquist F (2016). Total-evidence dating under the fossilized birth-death process. *Systematic Biology* 65: 228-249.

## Known limits

- The topology of living species is fixed within a run; date several trees to cover topological uncertainty.
- Fossils must fall inside a clade constraint below the root in root-conditioned runs.
- With rate shifts, the number of shifts is chosen by marginal likelihood across runs, not estimated within one run (no reversible jump). On large trees, estimated shift times can have several peaks; check them with two seeds, and use heated chains if the seeds disagree. There are no mass extinctions or extra sampling at the shifts. On small trees (30 to 120 living species), two shifts are rarely detected (test 2d).
- The result depends on the quality of the branch lengths and their standard errors.

## Citation and license

Not yet published. No license has been chosen yet, so please contact the author before reusing the code.
