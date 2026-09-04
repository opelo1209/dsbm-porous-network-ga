**Reviewer comment:**
"Implementation. Code via GitHub or similar platforms are required to cross validate results and implementations accuracy."

**Response:**
We thank the reviewer for this valuable suggestion. The source code has been made publicly
available on GitHub at **https://github.com/opelo1209/dsbm-porous-network-ga** to allow full
cross-validation of the results and implementation reported in the manuscript. The repository
provides:

- The complete constraint-preserving genetic algorithm as a single C source file implementing
  the population initialization, fitness evaluation, and crossover/mutation operators described
  in Sec. 3.1-3.5 (Fisher-Yates population initialization from a shared baseline network,
  per-category bond crossover, and per-category site-relocation mutation at a fixed 15%
  acceptance probability). The same source builds both the sequential baseline and the
  OpenMP-parallel version used throughout Sec. 4, selected at compile time via the `-fopenmp`
  flag, with no algorithmic differences between the two.
- A README documenting the exact compilation commands, compiler requirements, and the
  `-fopenmp` flag; the runtime and compile-time parameters corresponding to L, the overlap
  coefficient Omega, population size, thread count, and mutation probability (Sec. 4.1-4.3);
  and the random-number-generator seeding mechanism used by the implementation.
- Shell scripts and Python plotting scripts that reproduce the population-size sweep behind
  Table 1 (Sec. 4.2), the thread-scalability and lattice-size sweep behind Fig. 12 (Sec. 4.3),
  the sequential-vs-parallel comparison behind Fig. 14 (Sec. 4.4), and the overlap snapshots and
  spatial correlation function C(r) behind Fig. 15-18 (Sec. 4.5), together with the reference
  hardware/compiler configuration used to obtain the reported execution times (Sec. 4.1).
- An open-source MIT license permitting unrestricted reuse and independent verification by the
  reviewer or the broader community.

This repository is referenced in the revised manuscript and will remain publicly available upon
publication.
