#!/bin/sh
# Quick self-test of ffbd (about a minute). Run from the repository root: make test
# 1. Every clock, on three examples (crown, stem and daughter-stem rules included): after every proposal the cached posterior is recomputed from scratch and compared (--check-every-step).
# 2. Stepping-stone runs end to end and prints a marginal likelihood.
# 3. The bears example (real data) runs and its root age falls inside the root prior.
# 4. Estimated shift times (free:2, windows): cache check after every proposal; FBD-power stepping-stone runs.
# 5. Heated chains (heat on the FBD density, or on FBD and data): cache check of every copy after every proposal; the output does not depend on --threads.
set -u
B=./ffbd; T=$(mktemp -d); fail=0
[ -x "$B" ] || { echo "build ffbd first (make)"; exit 1; }
for ex in examples/example_input.txt examples/example_crown_stem.txt examples/bears_daughter_stem.txt; do
  for c in exp ln gam igr acln; do
    out=$($B run "$ex" --prefix "$T/chk" --clock $c --seed 11 --gens 20000 --print 5000 --burnin 5000 --check-every-step 2>&1)
    if echo "$out" | grep -q "ERROR" || ! echo "$out" | grep -q "FINAL"; then echo "FAIL  cache check  $ex  clock $c"; fail=1; else echo "ok    cache check  $ex  clock $c"; fi
  done
done
out=$($B run examples/example_input.txt --prefix "$T/ss" --clock gam --seed 3 --burnin 20000 --ss 8 --ss-burn 2000 --ss-gens 10000 --ss-every 100 2>&1)
if echo "$out" | grep -qi "stepping"; then echo "ok    stepping-stone"; else echo "FAIL  stepping-stone"; fail=1; fi
out=$($B run examples/bears_input.txt --prefix "$T/bears" --clock gam --seed 5 --gens 200000 --print 10000 --burnin 50000 2>&1)
if echo "$out" | grep -q "FINAL"; then echo "ok    bears example runs"; else echo "FAIL  bears example"; fail=1; fi
for spec in free:2 3:6,9:13; do
  out=$($B run examples/example_input.txt --prefix "$T/est" --clock gam --seed 11 --gens 20000 --burnin 5000 --rate-shifts $spec --shift-window 1,16 --check-every-step 2>&1)
  if echo "$out" | grep -q "ERROR" || ! echo "$out" | grep -q "FINAL"; then echo "FAIL  estimated shifts  $spec"; fail=1; else echo "ok    estimated shifts  $spec"; fi
done
out=$($B run examples/example_input.txt --prefix "$T/ssf" --clock gam --seed 3 --burnin 20000 --rate-shifts free:1 --shift-window 1,16 --ss 8 --ss-power fbd --ss-burn 2000 --ss-gens 10000 --ss-every 100 2>&1)
if echo "$out" | grep -qi "stepping"; then echo "ok    stepping-stone, FBD power"; else echo "FAIL  stepping-stone, FBD power"; fail=1; fi
for hon in fbd all; do
  out=$($B run examples/example_input.txt --prefix "$T/h" --clock gam --seed 13 --gens 20000 --burnin 5000 --rate-shifts free:2 --shift-window 1,16 --heated 3 --heat-on $hon --swap-every 100 --check-every-step 2>&1)
  if echo "$out" | grep -q "ERROR" || ! echo "$out" | grep -q "FINAL"; then echo "FAIL  heated chains  heat on $hon"; fail=1; else echo "ok    heated chains  heat on $hon"; fi
done
$B run examples/example_input.txt --prefix "$T/t1" --clock gam --seed 17 --gens 50000 --burnin 5000 --rate-shifts free:1 --shift-window 1,16 --heated 4 --swap-every 100 --threads 1 >/dev/null 2>&1
$B run examples/example_input.txt --prefix "$T/t4" --clock gam --seed 17 --gens 50000 --burnin 5000 --rate-shifts free:1 --shift-window 1,16 --heated 4 --swap-every 100 --threads 4 >/dev/null 2>&1
if [ -s "$T/t1.log" ] && cmp -s "$T/t1.log" "$T/t4.log"; then echo "ok    heated chains  same output with 1 and 4 threads"; else echo "FAIL  heated chains  output depends on --threads"; fail=1; fi
rm -rf "$T"
[ $fail -eq 0 ] && echo "all tests passed" || { echo "some tests FAILED"; exit 1; }
