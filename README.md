# Constraint-Preserving Genetic Algorithm for DSBM Porous Networks

Reference implementation for:

> Moreno-Montiel, B., Cordero-Sánchez, S., Villegas-Cortez, J. "A Constraint-Preserving
> Genetic Algorithm for the Construction of Porous Networks Based on the Dual Site-Bond Model."

This repository accompanies the paper for cross-validation of results and implementation
accuracy. It contains the genetic algorithm's C source (sequential and OpenMP-parallel, see
below), scripts to reproduce Table 1 and Figures 12-18, and a pure Monte Carlo baseline
implementation (Sec. 2.2's classical approach) used as a quantitative point of comparison.

**Please read "Reproducibility notes / known limitations" near the end of this file before
treating any script output as a literal match to the paper's reported numbers.** The public
interface (interactive `scanf` prompts, a compile-time population size, no fixed iteration
budget, no fixed random seed) differs in a few concrete ways from the internal experimental
harness used to produce the paper's figures, and this README is explicit about each gap and
the closest faithful reproduction available.

## Repository layout

```
src/genetico/     ConstructorRedes2D_C4_Genetico_Final.c + .vscode/  - the genetic algorithm
src/montecarlo/    ConstructorRedes2D_C4_MonteCarlo.c + .vscode/      - pure Monte Carlo baseline (Sec. 2.2)
scripts/           build.sh, env_info.sh, run_*.sh                   - build + experiment automation (bash)
analysis/          plot_*.py, requirements.txt                       - Python plotting for Table 1 / Fig. 12-18 / GA vs. MC
results/           (created by the scripts; not committed)           - raw logs, CSVs and generated figures
```

Each program lives in its own subfolder (with its own `.vscode/`) specifically so they can each be opened
directly as a VS Code workspace root and built/run/debugged with the "C/C++ Runner" extension's one-click
buttons without the two `main()` functions colliding in a single build. Open `src/genetico/` or
`src/montecarlo/` in VS Code (not the repository root) for that workflow; use `scripts/build.sh` /
the `run_*.sh` scripts below for the scripted reproducibility pipeline.

## Sequential and parallel code

The source contains `#pragma omp` directives (Sec. 3's genetic operators and the fitness
evaluation are parallelized). It is intentionally a **single file** for both configurations:

- Compiled **without** `-fopenmp`, the OpenMP pragmas are ignored by the compiler and the
  program runs single-threaded. This is the "sequential" baseline used in Fig. 14.
- Compiled **with** `-fopenmp`, the same pragmas are honored and the program runs
  multi-threaded, controlled at run time by the standard `OMP_NUM_THREADS` environment
  variable (or all logical cores if unset).

This matches how the paper itself describes the two versions ("the parallel code was designed
to maintain the same behavior ... as the sequential version").

## Requirements

- A C compiler with C99 and OpenMP support (GCC/MinGW-w64 or Clang). Tested informally with
  GCC via MinGW-w64 on Windows; any recent GCC/Clang on Linux/macOS should work identically.
- `bash` for the scripts (Git Bash or WSL on Windows, native on Linux/macOS).
- Python 3.9+ with the packages in `analysis/requirements.txt` for plotting.

Run `scripts/env_info.sh` and keep its output (compiler version, `-fopenmp` support check,
OS, CPU) alongside any results you share — see "Reproducibility notes" below.

```bash
scripts/env_info.sh > results/environment.txt
```

## Building

```bash
scripts/build.sh
```

produces `bin/sequential/ga_seq` (no `-fopenmp`) and `bin/parallel/ga_omp` (`-fopenmp`), both
built with `gcc -O3 -std=c99 -Wall`. Override the compiler with `CC=clang scripts/build.sh`.

Equivalent manual commands, if you prefer not to use the script:

```bash
gcc -O3 -std=c99 -Wall src/genetico/ConstructorRedes2D_C4_Genetico_Final.c -o ga_seq -lm
gcc -O3 -std=c99 -Wall -fopenmp src/genetico/ConstructorRedes2D_C4_Genetico_Final.c -o ga_omp -lm
```

## Pure Monte Carlo baseline

`src/montecarlo/ConstructorRedes2D_C4_MonteCarlo.c` implements the classical baseline the paper contrasts
the GA against (Sec. 1.2/2.2, citing Cruz et al.): lattice elements are repeatedly exchanged
through random site-site and bond-bond permutations until the Construction Principle is
satisfied, with **no** population, crossover, mutation, size-category restriction, or relaxation
fallback. A candidate exchange is accepted only if it does not increase the total error
(greedy/blind acceptance, not simulated annealing); it shares the exact same node layout, Type 1
error accounting (`calcularErrT1Local`), and CSV export format as the genetic algorithm, so
results are directly comparable and the existing `plot_site_distribution.py` /
`plot_correlation.py` scripts work on its output unchanged. Naively re-evaluating the whole
lattice after every single exchange would be O(L^2) per attempt and intractable at the paper's
lattice sizes; instead only the (at most four) positions whose error can change as a result of a
given exchange are recomputed, an O(1)-per-attempt algorithm documented in the source.

`scripts/build.sh` builds it to `bin/montecarlo/mc_puro` alongside the GA binaries. It is
interactive (`scanf`) and, unlike the GA binary, takes **no** command-line arguments and has
**no** `argv` quirk:

```bash
printf '50\n410\n400\n50\n20000000\n' | ./bin/montecarlo/mc_puro
```

in order: `L`, `mediaS`, `mediaE`, `desviacion`, and a maximum number of attempts (`0` = no
limit; use with caution - see below). Because there is no relaxation fallback, convergence is
**not guaranteed** within any fixed budget for large lattices - this is expected and is the
whole point of using it as the "blind search" baseline; the program reports whichever residual
error remains if the attempt budget is exhausted, rather than hanging indefinitely.

`scripts/run_ga_vs_mc_comparison.sh` runs both implementations on the same lattice sizes and
`(mediaS, mediaE, desviacion)` and produces a side-by-side execution-time / attempts-to-converge
comparison, plotted with `python analysis/plot_ga_vs_mc.py results/ga_vs_mc/ga_vs_mc_raw.csv`.

## Running

The GA program is interactive (`scanf`) and takes four inputs in this order:

1. Lattice size `L` (creates an `L x L` periodic lattice, Sec. 3.1)
2. Mean site radius `mediaS`
3. Mean bond radius `mediaE` (must be `< mediaS` for a physically meaningful overlap)
4. Half-width / standard deviation `desviacion` (site and bond radii are drawn uniformly
   from `[media - desviacion, media + desviacion]`)

**Known quirk:** the program unconditionally reads `argv[2]` (reserved for an unused,
never-implemented distribution-selection flag) before doing anything with the four `scanf`
inputs above. Always invoke it with at least two command-line arguments, even if they are
ignored placeholders:

```bash
./bin/parallel/ga_omp ignored ignored
# then type: L, mediaS, mediaE, desviacion (one per line, or piped as below)
```

Example non-interactive invocation (population = 64, Omega = 0.9 per Sec. 4.2/4.3):

```bash
export OMP_NUM_THREADS=32
printf '50\n410\n400\n50\n' | ./bin/parallel/ga_omp ignored ignored
```

The program prints its progress (`Optimo[gen] = fitness`) to stdout as new best solutions are
found, and on completion writes `red_colores.csv` (columns `x,y,r_Sitio,color`) to the current
working directory, plus a final `Tiempo total de ejecución: <seconds> segundos` line.

## Parameters used (L, Omega, population, iterations, threads, mutation probability)

| Parameter | Where it lives | Paper value(s) | Notes |
|---|---|---|---|
| `L` (lattice size) | `scanf` input #1 | 50 to 500 | Direct 1:1 runtime input. |
| Population size (`numCromosomas`) | **compile-time constant** in `main()`, `src/...c` | 4, 8, 16, 32, 64, 128 (64 optimal, Sec. 4.2) | Not a runtime input in this source. `scripts/run_population_sweep.sh` / `run_thread_sweep.sh` vary it by patching a *temporary* copy of the source with `sed` before compiling — `src/` itself is never modified. |
| Thread count | `OMP_NUM_THREADS` env var (parallel build only) | 4, 8, 16, 32, 64, 128 (32 optimal, Sec. 4.3) | Standard OpenMP mechanism; no source change needed. |
| Mutation probability | hard-coded `const double probMutacion = 0.15;` in `mutarCromosoma()` | 15% (Sec. 3.5) | Fixed in source, matches the paper exactly. |
| Iterations / generations | **fixed budget of 13000** via `MAX_GENERACIONES` in `main()` (0 = unbounded, run until `mejorFitness == 0`) | fixed budget of 13000 (Sec. 4.1/4.2) | Matches the paper by default. See "Reproducibility notes" item 1. |
| Omega (overlap) | **not a direct input** — implied by `mediaS`, `mediaE`, `desviacion` | 0.3, 0.6, 0.9 | See the derivation and concrete triples below. |
| Random seed | `srand(time(NULL))` in `main()` | not reported per-run in the paper | See "Random seeds" below. |

### Omega -> (mediaS, mediaE, desviacion)

The current interface samples both site and bond radii uniformly from a window of half-width
`desviacion` around their mean. Treating both as same-width uniform windows, the fraction of
one window covered by the other is:

```
Omega = 1 - |mediaS - mediaE| / (2 * desviacion)
```

This reproduces `Omega = 0.9` exactly for `(mediaS=410, mediaE=400, desviacion=50)`, the
values Sec. 4.2/4.3 of the paper state explicitly, and was cross-checked against two
independently logged runs of an earlier iteration of this codebase that printed Omega
directly: `(440, 400, 50) -> Omega = 0.6` and `(440, 370, 50) -> Omega = 0.3`. The three
concrete triples used by `scripts/run_overlap_snapshots.sh` are:

| Omega | mediaS | mediaE | desviacion |
|---|---|---|---|
| 0.3 | 440 | 370 | 50 |
| 0.6 | 440 | 400 | 50 |
| 0.9 | 410 | 400 | 50 |

### Random seeds

`main()` calls `srand(time(NULL))` once at startup and exposes no fixed-seed option, so:

- **Exact, bit-for-bit reproduction of a single run's output is not possible** with the
  current public interface.
- **Statistical reproducibility** (convergence trends, execution-time scaling, the qualitative
  effect of Omega on spatial correlation) is expected to hold across independent runs — this
  is consistent with how the paper itself reports results as mean ± SD over repeated trials
  (Table 1), rather than as a single fixed-seed run.
- If bit-exact reproduction is required, the one-line change needed is replacing
  `srand(time(NULL))` with `srand(SEED)` for a fixed integer `SEED` in `src/...c`. This is
  intentionally left undone here rather than silently changing the algorithm's shipped
  behavior; treat it as a documented extension point.

## Reproducing Table 1 and Figures 12-18

All experiment scripts live in `scripts/` and write raw CSVs/logs under `results/`; the
matching `analysis/plot_*.py` script turns each CSV into a figure. Install the Python
dependencies once:

```bash
python -m venv .venv && source .venv/bin/activate   # optional
pip install -r analysis/requirements.txt
```

| Paper artifact | Script | Plot |
|---|---|---|
| Table 1 (population-size sweep) | `scripts/run_population_sweep.sh` | `python analysis/plot_table1.py results/table1/table1_convergence.csv` |
| Fig. 12 (execution time vs. L, per thread count) | `scripts/run_thread_sweep.sh` | `python analysis/plot_execution_time.py results/fig12_13/fig12_13_raw.csv` |
| Fig. 13 (minimum error vs. L, per thread count) | `scripts/run_thread_sweep.sh` (same run) | `python analysis/plot_min_error.py results/fig12_13/fig12_13_raw.csv` (see caveat below — will be flat at 0) |
| Fig. 14 (sequential vs. parallel execution time) | `scripts/run_seq_vs_parallel.sh` | `python analysis/plot_seq_vs_parallel.py results/fig14/fig14_raw.csv` |
| Fig. 15-17 (site-size snapshots at Omega = 0.3/0.6/0.9) | `scripts/run_overlap_snapshots.sh` | `python analysis/plot_site_distribution.py <csv> --title "Omega = 0.X"` |
| Fig. 18 (spatial correlation C(r), log-log) | same run as Fig. 15-17 | `python analysis/plot_correlation.py <csv_0.3> <csv_0.6> <csv_0.9>` |
| GA vs. pure Monte Carlo baseline (Sec. 1.2/2.2, quantitative) | `scripts/run_ga_vs_mc_comparison.sh` | `python analysis/plot_ga_vs_mc.py results/ga_vs_mc/ga_vs_mc_raw.csv` |

Every `run_*.sh` script has a header comment documenting its exact parameters and how to
override them (population size, L values, thread counts, timeout). **Large-L, few-thread
configurations are genuinely slow** — the paper itself reports ~28,000 s for the sequential
L=500 case — so the scripts default to a smaller, fast smoke-test range and document the
environment variable to widen it to the paper's full range.

Quick smoke test (a few seconds, just to confirm the toolchain works end-to-end):

```bash
scripts/build.sh
printf '20\n410\n400\n50\n' | ./bin/parallel/ga_omp ignored ignored
python analysis/plot_site_distribution.py red_colores.csv --title "Smoke test"
```

## Reproducibility notes / known limitations

These are intentional properties of the current public source (left unmodified from what was
already implemented/corrected in this repository), documented here rather than silently
patched, so the gap between "what the paper reports" and "what this code will print if you
run it" is explicit:

1. **Fixed 13000-generation budget by default**, matching the paper (Sec. 4.1/4.2, Table 1/Fig.
   13). `main()`'s loop is `while (mejorFitness != 0 && (MAX_GENERACIONES <= 0 || gen <
   MAX_GENERACIONES))`, with `MAX_GENERACIONES = 13000`. If the budget runs out before
   `mejorFitness` reaches 0, the loop stops, the best network found so far is still exported to
   `red_colores.csv`, and the program prints the residual error instead of "converged" - matching
   how the paper itself reports Table 1/Fig. 13 (a fixed budget, generally non-zero residual
   error), rather than forcing exact convergence. This matters most for high-overlap
   configurations (e.g. Omega=0.9, `mediaS`/`mediaE` close together) where reaching error=0 via
   crossover/mutation alone can take an impractically long time even with population=64 and the
   guards below. Set `MAX_GENERACIONES = 0` in `main()` to restore the original unbounded
   behavior (run until `mejorFitness == 0`, however long that takes). Two additional mechanisms
   help close the distance to zero (or to whatever the budget allows) faster: a relaxation
   fallback (`relajaSitios`) once the residual violation rate drops below 0.1% (widened from an
   original 0.001% threshold that could never trigger below L~317, i.e. for every L the paper's
   own experiments use - see item 9 below), and a stagnation guard that refreshes part of the
   population if the best fitness stalls for too long (see item 10). `scripts/run_population_sweep.sh`
   and `plot_table1.py` currently report generations/time **to full convergence** rather than
   residual error at a fixed budget; with `MAX_GENERACIONES` now enforced, they could be updated
   to report residual error directly instead, closer to the paper's own Table 1 format.
2. **Population size is a compile-time constant** (`numCromosomas` in `main()`), not a runtime
   argument. The sweep scripts vary it by compiling a `sed`-patched temporary copy of the
   source per data point; `src/` itself is never modified.
3. **Omega is not a direct input or output.** It is implied by `(mediaS, mediaE, desviacion)`;
   see the derivation above.
4. **Random seed is time-based**, not fixed — see "Random seeds" above.
5. **`argv[2]` is read unconditionally** even though it is unused — always pass at least two
   command-line arguments (see "Running" above).
6. **Fitness accounting detail:** `evaluarFitness()` checks each of a site's four incident
   bonds against *that* site; summed over the whole lattice, every physical bond therefore gets
   evaluated once from **each** of its two endpoint sites, so a bond that violates both
   neighbors can contribute up to 2 to the total error, whereas Eq. 2 in the paper defines a
   single OR-based indicator per bond (contributing at most 1). Both formulations reach exactly
   0 under the same "fully consistent lattice" condition, so this does not affect the
   convergence criterion, only the intermediate magnitude of the reported error.
7. **Fixed:** `relajaSitios` (the legacy relaxation fallback, Sec. 2.2) used to loop over
   `numCromosomas` iterations of the *same* single chromosome it is given (its `NODO_BSM **`
   parameter did not vary with the loop index), applying its `+0.1%` site-radius growth
   `numCromosomas` times per call instead of once, and doing so under an OpenMP
   `collapse(3)` with no synchronization between the redundant iterations (a data race on
   `r_Sitio`). The redundant loop has been removed; the growth is now applied exactly once
   per call, deterministically, in both the sequential and parallel builds.

8. **The pure Monte Carlo baseline has no convergence guarantee within a bounded budget**, by
   design (see "Pure Monte Carlo baseline" above) - unlike the GA, it has no relaxation fallback,
   so `scripts/run_ga_vs_mc_comparison.sh` may legitimately report a non-zero residual error for
   it at large L once its `--max-attempts` cap is hit. This is not a bug; it is the comparison's
   point. Its output also tends to show large, smoothly-bounded spatial domains of same-category
   sites (visibly, and confirmed numerically: ~78% same-category neighbor pairs on a run vs. ~33%
   expected under spatial randomness) rather than a well-mixed arrangement. This is not a
   plotting artifact - it is an expected consequence of its "never accept a worse exchange"
   acceptance rule, which behaves like a zero-temperature local-search dynamics known to produce
   domain coarsening (the same class of phenomenon behind phase separation in zero-temperature
   Ising/Glauber dynamics or voter models), and is itself a concrete illustration of the paper's
   point (Sec. 1.2/2.2) that blind Monte Carlo search does not guarantee a physically realistic
   spatial arrangement the way the GA's category-preserving operators do.

9. **`relajaSitios`'s trigger threshold was widened from 0.001% to 0.1%** (see item 1). The
   original `< 0.001` percentage comparison required `error < 0.025` for L=50 - impossible for
   an integer error count - so the fallback could never fire for any L below ~317, silently
   defeating its documented purpose (finishing off the last few violations near convergence)
   for every L value the paper's own experiments use (50-500... except the very largest are
   still below 317). At 0.1% it becomes reachable (`error <= 2` at L=50, `<= 250` at L=500)
   while still only firing once the population is genuinely close to a feasible solution.

10. **Stagnation guard (new, additive - not one of the paper's Sec. 3 operators).** The
    elitist replacement in `cruzarCromosomas` (keep the 2 best of 2 parents + 2 offspring) can
    collapse the population into near-identical chromosomes after enough generations, leaving
    mutation's single-site swaps as the only remaining source of novelty - which stalls for a
    long time on violations that need several coordinated changes at once. If `mejorFitness`
    hasn't improved for 50 generations, `inyectarDiversidad` replaces the worst-performing
    quarter of the population (by that generation's fitness, never the current best) with fresh
    Fisher-Yates permutations of the same `M_base` used at startup - identical mechanism to
    `inicializarPoblacion`, so `F_S(R_S)`/`F_B(R_B)` stay exactly invariant; only which spatial
    arrangements are present in the population changes. Prints
    `Estancamiento detectado: N individuos reemplazados...` when it fires.

None of the above affect the correctness of the core constraint-preserving genetic operators
(population initialization, crossover, mutation — Sec. 3.2/3.4/3.5), which this source
implements as described in the paper (Fisher-Yates population initialization from a shared
baseline network, per-category bond crossover, per-category site-relocation mutation at a
fixed 15% acceptance probability).

## Citation

See `CITATION.cff`.

## License

MIT — see `LICENSE`.
