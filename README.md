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

## Output

Running the program prints step-by-step results to the console and also writes them to text files in the working directory, generated fresh on every run:

- `Lab_1.txt` — output from the small hand-built examples in `Main.cpp`.
- `Special_large_matrix.txt` — output from the dedicated large `SpecialMatrix` test case.
- `Computational_Experiment_Results.txt` — full results table from `ComputationalExperiment::run_experiment()`, including per-configuration average relative error and the dimension/range dependence analysis.

> **Note:** these `.txt` files are program output, regenerated on every run — they're listed here for reference but aren't meant to be tracked in version control (see `.gitignore`).

## What the experiment measures

For each combination of **matrix type** (tridiagonal / special), **conditioning** (well- / ill-conditioned), **dimension**, and **value range**, the experiment:

1. Generates a random matrix `A` satisfying (or deliberately violating) diagonal dominance.
2. Generates a random exact solution vector `x_exact`.
3. Computes `d = A * x_exact`.
4. Solves for `x_computed` using the Thomas algorithm.
5. Reports the **average relative error** between `x_computed` and `x_exact` over multiple trials.

This lets you see, empirically, how conditioning, problem size, and the magnitude of matrix entries affect the numerical stability of the Thomas algorithm.
