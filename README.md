# Tridiagonal Solver (Thomas Algorithm) — C++

A C++ implementation of the **Thomas algorithm** (tridiagonal matrix algorithm) for solving systems of linear algebraic equations (SLAE), extended to a custom **"special" matrix** type with an added dense last row. The project also includes a computational experiment that measures how solution accuracy is affected by matrix conditioning, system size, and value range.

## Background

A **tridiagonal system** `Ax = d` has a matrix with nonzero entries only on the main diagonal and its immediate neighbors:

```
| b1  c1               |   |x1|   |d1|
| a2  b2  c2           |   |x2|   |d2|
|     a3  b3  c3       | * |x3| = |d3|
|        ...  ...  ... |   |..|   |..|
|              an  bn  |   |xn|   |dn|
```

Instead of general Gaussian elimination (O(n³)), the **Thomas algorithm** solves this in **O(n)** time using a forward sweep (computing coefficients `L` and `M`) followed by a backward substitution.

This project extends that idea to a **"special" matrix**: a tridiagonal matrix where the **last row is fully dense** instead of only having 2 nonzero entries. This breaks the standard Thomas recurrence, so the algorithm is adapted to fold the extra coupling from the dense last row into the forward sweep before back-substitution.

## Project Structure

| File | Purpose |
|---|---|
| `Vector.h` / `Vector.cpp` | A lightweight custom vector class with manual memory management (deep copy, 1-based bounds-checked indexing, `norm()`, random fill, file/console printing). |
| `TridiagonalMatrix.h` | Standard tridiagonal matrix: storage as 3 diagonals, matrix-vector multiplication, and the classic Thomas algorithm solver. |
| `SpecialMatrix.h` | Tridiagonal matrix with an additional dense last row. Implements a modified Thomas algorithm that accounts for the extra row, plus diagnostics (compatibility checks, step-by-step equivalence verification, solution accuracy checks). |
| `ComputationalExperiment.h` / `.cpp` | Generates random well-conditioned and ill-conditioned test systems across a range of sizes and value ranges, solves them, and reports average relative error. |
| `Main.cpp` | Demonstrates all of the above on small hand-built examples, then runs the full computational experiment. |

## Features

- **1-based indexing** throughout, matching the mathematical notation used in the course methodology.
- **Diagonal dominance check** (`|b_i| ≥ |a_i| + |c_i|`) used as the well-conditioned criterion for generating test matrices; violating it is used to construct ill-conditioned systems.
- **Verbose mode** for the Thomas algorithm, printing every intermediate `L[i]` / `M[i]` coefficient to both console and file — useful for tracing the algorithm step by step.
- **Diagnostics for `SpecialMatrix`**:
  - `check_compatibility_conditions` — validates the last-row constraints.
  - `check_equivalence_step_by_step` — confirms the transformed system is equivalent to the original.
  - `check_solution_accuracy_verbose` — computes and prints the residual `‖Ax − d‖`.
- **Computational experiment** across dimensions `4 … 4096` and value ranges `1e-1 … 1e3`, comparing:
  - well-conditioned vs. ill-conditioned tridiagonal systems,
  - well-conditioned vs. ill-conditioned special (dense-last-row) systems,
  - how average relative error scales with dimension (fixed range) and with value range (fixed dimension).

## Build & Run

The project is set up as a Visual Studio solution (`lab_1.sln`, `lab_1.vcxproj`).

1. Open `lab_1.sln` in **Visual Studio** (2019+ recommended, C++17 or later).
2. Select **Debug** or **Release**, platform **x64**.
3. Build and run (`Ctrl+F5`).

No external dependencies are required — everything uses the C++ standard library (`<vector>`, `<random>`, `<fstream>`, etc.).

## Output files

Running the program prints step-by-step results to the console and also writes them to text files in the working directory, generated fresh on every run:

- `Lab_1.txt` — output from the small hand-built examples in `Main.cpp`.
- `Special_large_matrix.txt` — output from the dedicated large `SpecialMatrix` test case.
- `Computational_Experiment_Results.txt` — full results table from `ComputationalExperiment::run_experiment()`, including per-configuration average relative error and the dimension/range dependence analysis.

> **Note:** these `.txt` files are program output, regenerated on every run — they aren't tracked in version control (see `.gitignore`).

## Sample Results

### Correctness check (small tridiagonal & special matrices)

Both solvers are verified on a small hand-built 4×4/5×5 example, computing `X` from `d = A·X_exact` and comparing back:

```
Expected solution:      [1, -1, 2, -1, 1]
X_computed:             [1, -1, 2, -1, 1]
Max error:              0
Residual (A*x - d):     [0, 0, 0, 0, 0]
Residual norm (Ax-d):   0
```

### Large special matrix test (size 12)

For a randomly generated, well-conditioned 12×12 special matrix, the modified Thomas algorithm reproduces the exact solution to within floating-point precision:

```
A*x =        [39.339, 10.0791, -37.837, 53.1398, -55.5455, -15.1022, 62.6149, -27.5037, 52.8715, 29.183, -16.2757, -91.2133]
d   =        [39.339, 10.0791, -37.837, 53.1398, -55.5455, -15.1022, 62.6149, -27.5037, 52.8715, 29.183, -16.2757, -91.2133]
Residual A*x - d:  ~1e-14 (elementwise)
Residual norm:      1.42109e-14
error_norm:         6.11511e-13
```

### Computational experiment — average relative error

**Well-conditioned systems** (both tridiagonal and special matrices) stay near machine precision regardless of size or value range:

| Dimension | Range | Tridiagonal (well) | Special (well) |
|---|---|---|---|
| 4 | 1e-1 | 1.32e-16 | 3.22e-16 |
| 64 | 1e1 | 4.59e-16 | 1.03e-15 |
| 512 | 1e2 | 2.95e-16 | 5.52e-16 |
| 4096 | 1e3 | 3.98e-16 | 3.77e-15 |

**Ill-conditioned systems** show error growing by several orders of magnitude as dimension increases, especially at larger value ranges:

| Dimension | Range | Tridiagonal (ill) | Special (ill) |
|---|---|---|---|
| 4 | 1e-1 | 1.48e-15 | 9.96e-15 |
| 64 | 1e1 | 1.02e-13 | 4.74e-14 |
| 512 | 1e2 | 1.86e-11 | 1.16e-12 |
| 4096 | 1e3 | 2.97e-10 | 3.06e-06 |

### Dimension & range dependence (fixed special matrix, range = 1.0)

| Dimension | Well-conditioned error | Ill-conditioned error |
|---|---|---|
| 4 | 4.69e-15 | 1.69e-14 |
| 64 | 1.12e-15 | 6.50e-12 |
| 1024 | 1.19e-15 | 5.75e-11 |
| 4096 | 6.31e-16 | 3.76e-07 |

| Range | Well-conditioned error | Ill-conditioned error |
|---|---|---|
| 1e-1 | 7.63e-16 | 2.25e-12 |
| 1e1 | 1.68e-15 | 2.66e-12 |
| 1e2 | 5.10e-16 | 1.15e-13 |
| 1e3 | 8.03e-16 | 1.81e-12 |

### Overall average errors

```
Average error of well-conditioned tridiagonal matrix:  3.716641e-16
Average error of ill-conditioned tridiagonal matrix:   3.315236e-10
Average error of well-conditioned special matrix:      2.459164e-15
Average error of ill-conditioned special matrix:       8.171599e-08
```

### Key takeaway

Well-conditioned systems solve to **machine precision** regardless of matrix size or value range — the Thomas algorithm is numerically stable for diagonally dominant systems. Ill-conditioned systems, by contrast, accumulate error that grows by **several orders of magnitude** as the system size increases (roughly 1e-15 → 1e-7 across the tested range of dimensions), showing the expected numerical sensitivity of the method once diagonal dominance is violated.
