# FEM from Scratch

[![FEM Analysis](https://github.com/Vaib2001/fem-from-scratch/actions/workflows/fem-analysis.yml/badge.svg)](https://github.com/Vaib2001/fem-from-scratch/actions/workflows/fem-analysis.yml)

Finite Element Method implementations developed from first principles in Python and verified against analytical solutions.

The project focuses on the numerical mechanics behind FEM: element formulation, global assembly, boundary conditions, solution procedures, stress recovery, verification, mesh convergence, and automated numerical testing.

---

## Overview

The repository currently contains two one-dimensional structural mechanics problems.

The **uniform axial bar** provides a clean introduction to the FEM formulation and demonstrates analytical verification.

The **tapered axial bar** introduces a spatially varying cross-section, producing a nonlinear displacement field and a genuine finite-element discretization error. This makes it possible to investigate mesh convergence and estimate the numerical convergence rate.

### Current capabilities

- 2-node linear axial bar elements
- global stiffness-matrix assembly
- essential and natural boundary conditions
- nodal displacement solution
- strain and stress recovery
- analytical verification
- variable cross-sectional area
- mesh-refinement studies
- convergence-order estimation
- engineering visualization
- automated numerical tests
- continuous integration with GitHub Actions

---

# 1. Uniform Axial Bar

Consider a prismatic elastic bar with length \(L\), cross-sectional area \(A\), and Young's modulus \(E\).

The left end is fixed and an axial load \(P\) is applied at the right end.

```text
x = 0                                      x = L
│                                             →
│============================================= P
│
Fixed
```

The governing equation is

$$
-\frac{d}{dx}
\left(
EA\frac{du}{dx}
\right)
=0
$$

with boundary conditions

$$
u(0)=0
$$

and

$$
EA\frac{du}{dx}\bigg|_{x=L}=P.
$$

---

## Finite-element formulation

For a two-node linear bar element,

$$
\mathbf{u}^{(e)}
=
\begin{bmatrix}
u_1 \\
u_2
\end{bmatrix}.
$$

Using linear shape functions, the element stiffness matrix becomes

$$
\mathbf{k}^{(e)}
=
\frac{EA}{L_e}
\begin{bmatrix}
1 & -1 \\
-1 & 1
\end{bmatrix}.
$$

The element contributions are assembled into the global system

$$
\mathbf{K}\mathbf{u}
=
\mathbf{F}.
$$

After imposing the fixed displacement boundary condition, the reduced system is solved using `numpy.linalg.solve`.

---

## Analytical solution

For the uniform bar,

$$
u(x)
=
\frac{Px}{EA}.
$$

The axial strain is

$$
\varepsilon
=
\frac{du}{dx}
=
\frac{P}{EA},
$$

and therefore

$$
\sigma
=
E\varepsilon
=
\frac{P}{A}.
$$

The reference problem uses:

| Parameter | Value |
|---|---:|
| Length \(L\) | 1.0 m |
| Cross-sectional area \(A\) | 0.01 m² |
| Young's modulus \(E\) | 210 GPa |
| Applied load \(P\) | 100 kN |
| Finite elements | 4 |

The analytical tip displacement is

$$
u(L)
=
4.7619\times10^{-5}\ \text{m}
=
0.047619\ \text{mm},
$$

while the axial stress is

$$
\sigma
=
10\ \text{MPa}.
$$

---

## Displacement verification

![Uniform bar displacement](results/displacement_comparison.png)

The FEM nodes coincide with the analytical displacement field.

This is expected: the exact solution is linear in \(x\), while the two-node bar element also uses linear displacement interpolation. The exact solution therefore lies within the finite-element approximation space.

---

## Stress verification

![Uniform bar stress](results/stress_distribution.png)

The finite-element stress remains constant along the bar and agrees with the analytical value

$$
\sigma = 10\ \text{MPa}.
$$

---

# 2. Tapered Axial Bar

The second problem introduces a cross-sectional area that varies continuously along the bar.

The area is defined as

$$
A(x)
=
A_0
\left(
1+\alpha\frac{x}{L}
\right),
$$

where \(A_0\) is the area at the fixed end and \(\alpha\) controls the taper.

For the implemented example,

$$
\alpha=1,
$$

which gives

$$
A(L)=2A_0.
$$

Unlike the uniform-bar problem, the exact displacement field is no longer linear.

---

## Analytical solution

For a constant axial force \(P\),

$$
\frac{du}{dx}
=
\frac{P}{EA(x)}.
$$

Substituting the linearly varying area gives

$$
\frac{du}{dx}
=
\frac{P}
{
EA_0
\left(
1+\alpha x/L
\right)
}.
$$

Integrating from the fixed end gives

$$
u(x)
=
\frac{PL}
{
EA_0\alpha
}
\ln
\left(
1+\alpha\frac{x}{L}
\right).
$$

The exact stress field is

$$
\sigma(x)
=
\frac{P}{A(x)}.
$$

Because the exact displacement is logarithmic rather than linear, linear finite elements only approximate the solution.

---

## Tapered-bar displacement

![Tapered bar displacement](results/tapered_displacement.png)

The finite-element solution follows the analytical displacement field closely while retaining the piecewise-linear nature of the numerical approximation.

---

## Tapered-bar stress

![Tapered bar stress](results/tapered_stress.png)

As the cross-sectional area increases toward the free end, the axial stress decreases according to

$$
\sigma(x)
=
\frac{P}{A(x)}.
$$

The FEM element stresses approach the continuous analytical stress field.

---

# Mesh Convergence

A central requirement of a reliable finite-element implementation is that the numerical solution approaches the exact solution as the mesh is refined.

The tapered bar is solved using

$$
N =
1,\ 2,\ 4,\ 8,\ 16,\ 32,\ 64
$$

elements.

The relative tip-displacement error is evaluated as

$$
e_{\mathrm{rel}}
=
\frac{
\left|
u_{\mathrm{FEM}}(L)
-
u_{\mathrm{exact}}(L)
\right|
}{
\left|
u_{\mathrm{exact}}(L)
\right|
}.
$$

![Mesh convergence](results/tapered_convergence.png)

The error decreases systematically as the number of elements increases.

The observed behaviour is approximately

$$
e
\propto
N^{-2},
$$

which is consistent with the expected convergence behaviour of linear finite elements for this displacement quantity.

This provides a numerical verification that the implementation converges toward the analytical solution under mesh refinement.

---

# Automated Verification

The repository contains automated numerical tests using `pytest`.

The test suite verifies the fixed boundary condition, the analytical tip displacement of the uniform bar, the analytical uniform stress state, reduction of tapered-bar error under mesh refinement, and accuracy of a sufficiently refined tapered-bar solution.

The complete test suite is executed automatically by GitHub Actions after changes are pushed to the repository.

A successful workflow therefore checks both that the implementation executes correctly and that key numerical and physical properties remain satisfied.

---

# Continuous Integration

The GitHub Actions workflow automatically performs the complete analysis pipeline:

```text
Push to main
     │
     ▼
Create clean Python environment
     │
     ▼
Install dependencies
     │
     ▼
Run uniform-bar solver
     │
     ▼
Generate uniform-bar figures
     │
     ▼
Run tapered-bar solver
     │
     ▼
Run mesh-convergence study
     │
     ▼
Run numerical verification tests
     │
     ▼
Update generated FEM results
```

This makes the numerical analysis reproducible independently of the local development environment.

---

# Project Structure

```text
fem-from-scratch/
│
├── .github/
│   └── workflows/
│       └── fem-analysis.yml
│
├── src/
│   ├── __init__.py
│   ├── bar_1d.py
│   ├── tapered_bar.py
│   ├── visualize_bar.py
│   └── convergence_study.py
│
├── tests/
│   └── test_fem.py
│
├── results/
│   ├── displacement_comparison.png
│   ├── stress_distribution.png
│   ├── tapered_displacement.png
│   ├── tapered_stress.png
│   └── tapered_convergence.png
│
├── requirements.txt
├── pytest.ini
├── README.md
└── LICENSE
```

---

# Running the Project

Install the required packages:

```bash
pip install -r requirements.txt
```

Run the uniform-bar analysis:

```bash
python src/bar_1d.py
```

Generate the uniform-bar figures:

```bash
python src/visualize_bar.py
```

Run the tapered-bar analysis:

```bash
python src/tapered_bar.py
```

Run the convergence study:

```bash
python src/convergence_study.py
```

Run the numerical verification tests:

```bash
pytest -q
```

---

# Engineering Concepts Demonstrated

This repository demonstrates the progression from an elementary finite-element formulation to numerical verification of a non-trivial problem.

In particular, it covers stiffness formulation, matrix assembly, degree-of-freedom handling, displacement boundary conditions, force boundary conditions, stress recovery, analytical verification, discretization error, mesh refinement, convergence-rate estimation, scientific visualization, automated testing, and continuous integration.

---

# Future Development

The next stage of the project will extend the formulation from one-dimensional axial mechanics toward multi-dimensional structural FEM.

Planned extensions include spring elements, 2D truss elements, local-to-global coordinate transformations, arbitrary truss geometries, distributed loading, 2D linear elasticity, triangular elements, and stress-field visualization.

---

## Technology

`Python` · `NumPy` · `Matplotlib` · `pytest` · `GitHub Actions`

---

## License

This project is released under the MIT License.
