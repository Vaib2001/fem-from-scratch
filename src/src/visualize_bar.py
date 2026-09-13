from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bar_1d import (
    analytical_displacement,
    analytical_stress,
    solve_bar,
)


def main():
    # ---------------------------------------------------------
    # Engineering problem
    # ---------------------------------------------------------

    length = 1.0
    area = 0.01
    youngs_modulus = 210e9
    force = 100000.0
    num_elements = 4

    # ---------------------------------------------------------
    # FEM solution
    # ---------------------------------------------------------

    nodes, displacements, stresses = solve_bar(
        length,
        area,
        youngs_modulus,
        force,
        num_elements,
    )

    # ---------------------------------------------------------
    # Analytical solution
    # ---------------------------------------------------------

    x_exact = np.linspace(
        0.0,
        length,
        300,
    )

    u_exact = analytical_displacement(
        x_exact,
        area,
        youngs_modulus,
        force,
    )

    sigma_exact = analytical_stress(
        area,
        force,
    )

    # ---------------------------------------------------------
    # Results directory
    # ---------------------------------------------------------

    project_root = Path(__file__).resolve().parents[1]

    results_dir = project_root / "results"

    results_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ---------------------------------------------------------
    # Plot 1:
    # FEM vs analytical displacement
    # ---------------------------------------------------------

    plt.figure(
        figsize=(8, 5)
    )

    plt.plot(
        x_exact,
        u_exact * 1000,
        label="Analytical solution",
    )

    plt.plot(
        nodes,
        displacements * 1000,
        "o",
        label="FEM nodes",
    )

    plt.xlabel(
        "Position along bar [m]"
    )

    plt.ylabel(
        "Axial displacement [mm]"
    )

    plt.title(
        "1D Axial Bar: FEM vs Analytical Solution"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        results_dir
        / "displacement_comparison.png",
        dpi=300,
    )

    plt.close()

    # ---------------------------------------------------------
    # Plot 2:
    # Stress distribution
    # ---------------------------------------------------------

    element_centres = (
        nodes[:-1]
        + nodes[1:]
    ) / 2

    plt.figure(
        figsize=(8, 5)
    )

    plt.plot(
        element_centres,
        stresses / 1e6,
        "o-",
        label="FEM stress",
    )

    plt.axhline(
        sigma_exact / 1e6,
        linestyle="--",
        label="Analytical stress",
    )

    plt.xlabel(
        "Position along bar [m]"
    )

    plt.ylabel(
        "Axial stress [MPa]"
    )

    plt.title(
        "1D Axial Bar: Stress Distribution"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        results_dir
        / "stress_distribution.png",
        dpi=300,
    )

    plt.close()

    print(
        "Results saved to:",
        results_dir,
    )


if __name__ == "__main__":
    main()
