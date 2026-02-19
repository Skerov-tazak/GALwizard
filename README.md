# GALwizard 
[![Test](https://github.com/Skerov-tazak/GALwizard/actions/workflows/test.yml/badge.svg)](https://github.com/Skerov-tazak/GALwizard/actions/workflows/test.yml)

**GALwizard** is a numerical linear algebra and complex math library written from scratch in C++. 

I built this over the summer of 2023 and 2024 because I wanted to understand exactly what happens under the hood when we blindly call `numpy.linalg.solve()` or use BLAS or LAPACK, I implemented the decompositions and solvers by myself but chose to retain `std::vector` as I wanted to focus on the linear algebra aspects.

This project was a deep dive into numerical stability, algorithmic complexity, and writing C++ code that can be extended, revised and maintained without needing costly refactoring and rewriting the code base. I learnt a lot about operators, overloading, input/output streams, memory management, writing clean code and how much impact simple design decisions made at the very beginning can have on later performance when a new unexpected feature is added. 

## What's Inside?

### Matrix Math & Decompositions
* **Arithmetic:** Matrix addition, multiplication and hermitian form. 
* **QR Decomposition:** Implemented two ways—using **Gram-Schmidt Orthogonalization** and **Householder Reflections** (because numerical stability matters).
* **PLU Decomposition:** Extracts the P, L, and U matrices with partial pivoting to handle singular or near-singular matrices.
* **Givens Rotations:** Uses similarity transformations to reduce matrices into Upper Hessenberg form.

### Linear System Solvers
* **Ax = b systems:** This goes a step further by calculating both the **particular solution** and the **nullspace** (basis vectors) for systems with infinite solutions.
* **Projections:** Handles rank-deficient matrices and finds the best-fit projected solutions for overdetermined systems.
* **Gaussian tools:** Row Echelon and Reduced Row Echelon (RREF) forms.
* **Space Equality:** Compares the column spaces of two matrices.

### Custom Complex Number Engine
I wrote a custom `number` class from the ground up to handle complex arithmetic accurately:
* Native Cartesian/Polar (Euler) conversions.
* Handles floating-point drift using a custom `approx_equal` tolerance (default epsilon = 1e-6).
* Calculates complex exponentiation, logarithms, and n-th complex roots.

## Installation and Setup 

GALwizard uses CMake and requires a C++11 compatible compiler. The testing framework, Catch2, is fetched automatically by CMake during the build process..

### 1. Clone the repository
```bash
git clone https://github.com/Skerov-tazak/GALwizard.git
cd GALwizard
```
### 2. Build using CMake
```bash
cmake -B build -S . 
cmake --build build
```
Now you should now be able to find the shared library file and 3 executable files in the build folder. 
 - `gal_tests`   Runs Catch2 tests on the library.
 - `num_tests` Runs Catch2 tests on the complex number class and utilities.
 - `manual_tests` Executes code from `main.cpp`.
 - `libGALwizard.so` Functional shared library code.

### 3. Editor Setup
You probably use an LSP like `clangd` (via Coc, or VSCode IntelliSense or whatever). CMake is configured to automatically generate a `compile_commands.json` file.

To ensure your editor resolves the `#include` paths correctly without throwing "file not found" errors, just symlink the database to your project root:
```bash
ln -s build/compile_commands.json .
```

## Quick Start Example

Here is a quick look at how the library API works for a PLU Decomposition:

```cpp
#include "include/functions.h"
#include <iostream>

using namespace gal;

int main() {
    // Initialize a 3x3 matrix (zeros by default)
    Matrix A = Matrix::nullmatrix(3, 3);
    
    A[0][0] = 0; A[0][1] = 2; A[0][2] = 1;
    A[1][0] = 1; A[1][1] = 0; A[1][2] = 2;
    A[2][0] = 2; A[2][1] = 1; A[2][2] = 0;

    // Perform PLU Decomposition (handles the necessary row swaps)
    std::vector<Matrix> plu = PLU_decompose(A);
    Matrix P = plu[0];
    Matrix L = plu[1];
    Matrix U = plu[2];

    std::cout << "Upper Triangular Matrix (U): \n" << U << std::endl;
    
    return 0;
}
```

## Design Notes

While making this project, I made a few specific architectural choices. Some of them were thought through from the start, others could be easily implemented thanks to modular design - yet other features were really painful to introduce later, but were still added due to great performance payoff. 

-   **One Dimensional Layout:** Instead of a fragmented `vector<vector<number>>`, matrices are implemented by a single contiguous  `std::vector<number>`. Row and column operations compute indexes using a `real_index()` helper function. This guarantees spatial locality and cache-friendly memory access.

- **Column Major Order** I used (or came up with) what is called column major order - that is the one dimensional array comprises of fixed-length columns stitched together. This idea came from an observation that column appending, swapping and extracting operations are much more common than row ones. Due to a modular design was easy to implement - as this was not part of the project from its onset. 
    
-   **Proxy Objects for Intuitive Syntax:** To maintain the one dimensional array while allowing standard syntax like `A[row][col]`, the overloaded `operator[]` returns lightweight `RowProxy` and `ConstRowProxy` objects, which handle the operations seamlessly.
    
-   **Tolerance-Based Zeroing:** Floating-point arithmetic in linear algebra is notoriously messy. The custom `number` class implements an `approx_equal` function with a strict epsilon threshold (1e-6). The library actively clears small values to prevent floating-point errors from impacting algorithms.

## Roadmap

There are many ways in which the library could be improved. Some of the planned short-term features are:

- [ ] More Matrix decompositions (SVD, SDS, EVD)
- [ ] The solution to the Eigenvalue problem with QR iteration (Power Method Almost Finished)
- [ ] Least Squares support via QR decomposition 

And in the long term I believe it would be amazing to implement:

- [ ] **Matrix Inheritance and Polymorphism:** There are many types of matrices which may need special handling and care - like tridiagonal, hessenberg, sparse or positive definite matrices. These appear often in physical systems (With many tensors being positive definite for example). Because of that this is one of my greatest plans for this library.
- [ ] **Numerical Solutions to Polynomial Problems:** Many interesting and useful questions involving polynomials can be turned into linear algebra questions (Like Gradient Descent). Adding a polynomial system and implementing algorithms that operate on polynomials would increase the library's usefulness beyond linear algebra (Though polynomials are a linear space so maybe we aren't really venturing too far...).   
- [ ] **Parallelization:** Many if not most problems in Linear Algebra can be easily parallelized. I want to integrate parallelize matrix multiplications and introduce vectorization - as this would significantly increase performance.
- [ ] **Pybind:** Most useful scripts start as simple ideas that are then quickly written in python. I believe it is necessary to know how to transform C++ code into python functions. 

