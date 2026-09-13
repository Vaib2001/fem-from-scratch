# FEM from Scratch

[![FEM Analysis](https://github.com/Vaib2001/fem-from-scratch/actions/workflows/fem-analysis.yml/badge.svg)](https://github.com/Vaib2001/fem-from-scratch/actions/workflows/fem-analysis.yml)

Finite Element Method implementations developed from first principles in Python, with analytical verification and engineering-focused post-processing.

The goal of this repository is to demonstrate the numerical mechanics behind FEM rather than relying on a commercial finite-element solver.

---

## Current implementation

### 1D Axial Bar

The first model considers a uniform elastic bar:

- fixed at the left end,
- subjected to an axial point load at the right end,
- discretized using 2-node linear bar elements.

The implementation includes:

- mesh generation,
- element stiffness formulation,
- global stiffness assembly,
- boundary-condition application,
- solution of the reduced system,
- strain and stress recovery,
- analytical verification,
- numerical error calculation,
- displacement visualization,
- stress visualization.

---

## Governing equation

For a one-dimensional elastic bar,

$$
-\frac{d}{dx}
\left(
EA\frac{du}{dx}
\right)=0
$$

where:

- \(E\) = Young's modulus,
- \(A\) = cross-sectional area,
- \(u(x)\) = axial displacement.

The boundary conditions are

$$
u(0)=0
$$

and

$$
EA\frac{du}{dx}\bigg|_{x=L}=P.
$$

---

## Finite-element formulation

For a two-node linear bar element, the element stiffness matrix is

$$
\mathbf{k}^{(e)}
=
\frac{EA}{L_e}
\begin{bmatrix}
1 & -1 \\
-1 & 1
\end{bmatrix}.
$$

The individual element matrices are assembled into the global system

$$
\mathbf{K}\mathbf{u}=\mathbf{F}.
$$

After applying the prescribed displacement at the fixed boundary, the reduced system is solved using NumPy.

---

## Analytical solution

For a uniform bar subjected to an axial end load,

$$
u(x)=\frac{Px}{EA}.
$$

The axial stress is constant:

$$
\sigma=\frac{P}{A}.
$$

For the example problem,

| Parameter | Value |
|---|---:|
| Length \(L\) | 1.0 m |
| Area \(A\) | 0.01 m² |
| Young's modulus \(E\) | 210 GPa |
| Applied load \(P\) | 100 kN |
| Elements | 4 |

The analytical tip displacement is

$$
u(L)=4.7619\times10^{-5}\ \text{m}
$$

or approximately

$$
0.04762\ \text{mm}.
$$

The corresponding axial stress is

$$
\sigma=10\ \text{MPa}.
$$

---

## FEM vs analytical displacement

The FEM nodal solution lies directly on the analytical solution.

![FEM vs analytical displacement](results/displacement_comparison.png)

For this particular problem, the analytical displacement field is linear. Since the finite elements also use linear interpolation, even a coarse linear-element mesh can reproduce the displacement field essentially exactly.

---

## Stress distribution

The axial stress remains constant along the uniform bar and agrees with the analytical value of 10 MPa.

![Axial stress distribution](results/stress_distribution.png)

---

## Project structure

```text
fem-from-scratch/
├── .github/
│   └── workflows/
│       └── fem-analysis.yml
│
├── src/
│   ├── bar_1d.py
│   └── visualize_bar.py
│
├── results/
│   ├── README.md
│   ├── displacement_comparison.png
│   └── stress_distribution.png
│
├── requirements.txt
├── README.md
└── LICENSE
```

---

## Run locally

Clone the repository and install the dependencies:

```bash
pip install -r requirements.txt
```

Run the FEM analysis:

```bash
python src/bar_1d.py
```

Generate the result figures:

```bash
python src/visualize_bar.py
```

The generated figures are stored in the `results/` directory.

---

## Automated analysis

A GitHub Actions workflow automatically:

1. creates a clean Python environment,
2. installs the project dependencies,
3. executes the FEM solver,
4. generates the engineering plots,
5. stores the generated results as workflow artifacts.

This helps verify that the implementation remains reproducible outside the development machine.

---

## Skills demonstrated

This project currently demonstrates:

- Finite Element Method fundamentals
- element stiffness formulation
- global matrix assembly
- boundary-condition treatment
- linear-system solution
- strain and stress recovery
- analytical verification
- scientific computing with NumPy
- engineering visualization with Matplotlib
- reproducible workflows with GitHub Actions

---

## Planned extensions

The repository will progressively move from simple 1D elements toward more realistic structural mechanics problems:

- variable-area / tapered axial bar
- mesh-convergence analysis
- spring elements
- 2D truss elements
- global coordinate transformations
- distributed loading
- 2D elasticity
- triangular finite elements
- stress-field visualization

The next implementation will introduce a non-uniform bar so that mesh refinement produces a genuine convergence study.
