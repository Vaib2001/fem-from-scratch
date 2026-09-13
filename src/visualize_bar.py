from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bar_1d import (
    analytical_displacement,
    analytical_stress,
    solve_bar,
)


def apply_github_dark_theme():
    """
    Apply a GitHub-dark inspired plotting theme.
    """

    plt.rcParams.update({
        "figure.facecolor": "#0d1117",
        "axes.facecolor": "#161b22",
        "axes.edgecolor": "#30363d",
        "axes.labelcolor": "#c9d1d9",
        "axes.titlecolor": "#f0f6fc",
        "xtick.color": "#c9d1d9",
        "ytick.color": "#c9d1d9",
        "text.color": "#c9d1d9",
        "grid.color": "#30363d",
        "grid.linestyle": "--",
        "grid.linewidth": 0.8,
        "legend.facecolor": "#161b22",
        "legend.edgecolor": "#30363d",
        "legend.framealpha": 1.0,
        "font.size": 12,
        "axes.titleweight": "bold",
        "axes.labelweight": "medium",
    })


def main():
    apply_github_dark_theme()

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

    x_exact = np.linspace(0.0, length, 300)

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

    tip_exact = analytical_displacement(
        length,
        area,
        youngs_modulus,
        force,
    )

    relative_error = abs(
        displacements[-1] - tip_exact
    ) / abs(tip_exact)

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

    fig, ax = plt.subplots(figsize=(9, 5.5))

    ax.plot(
        x_exact,
        u_exact * 1000,
        color="#58a6ff",
        linewidth=2.8,
        label="Analytical solution",
    )

    ax.plot(
        nodes,
        displacements * 1000,
        linestyle="None",
        marker="o",
        markersize=8,
        markerfacecolor="#f0883e",
        markeredgecolor="#f0f6fc",
        markeredgewidth=1.0,
        label="FEM nodes",
    )

    ax.set_title("1D Axial Bar: FEM vs Analytical Displacement", pad=14)
    ax.set_xlabel("Position along bar [m]")
    ax.set_ylabel("Axial displacement [mm]")

    ax.grid(True, alpha=0.8)
    ax.legend(loc="upper left")

    validation_text = (
        f"Tip displacement = {displacements[-1]*1000:.5f} mm\n"
        f"Analytical = {tip_exact*1000:.5f} mm\n"
        f"Relative error = {relative_error:.2e}"
    )

    ax.text(
        0.98,
        0.06,
        validation_text,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=10.5,
        bbox=dict(
            boxstyle="round,pad=0.4",
            facecolor="#0d1117",
            edgecolor="#30363d",
        ),
    )

    plt.tight_layout()

    plt.savefig(
        results_dir / "displacement_comparison.png",
        dpi=300,
        bbox_inches="tight",
        facecolor=fig.get_facecolor(),
    )

    plt.close(fig)

    # ---------------------------------------------------------
    # Plot 2:
    # Stress distribution
    # ---------------------------------------------------------

    element_centres = (
        nodes[:-1] + nodes[1:]
    ) / 2

    fig, ax = plt.subplots(figsize=(9, 5.5))

    ax.plot(
        element_centres,
        stresses / 1e6,
        color="#3fb950",
        linewidth=2.4,
        marker="o",
        markersize=8,
        markerfacecolor="#3fb950",
        markeredgecolor="#f0f6fc",
        markeredgewidth=1.0,
        label="FEM stress",
    )

    ax.axhline(
        sigma_exact / 1e6,
        color="#58a6ff",
        linestyle="--",
        linewidth=2.2,
        label="Analytical stress",
    )

    ax.set_title("1D Axial Bar: Axial Stress Distribution", pad=14)
    ax.set_xlabel("Position along bar [m]")
    ax.set_ylabel("Axial stress [MPa]")

    # Better limits so the constant stress line does not look cramped
    sigma_mpa = sigma_exact / 1e6
    ax.set_ylim(
        sigma_mpa - 0.6,
        sigma_mpa + 0.6,
    )

    ax.grid(True, alpha=0.8)
    ax.legend(loc="upper right")

    stress_text = (
        f"Stress = {sigma_mpa:.2f} MPa\n"
        f"Uniform bar\n"
        f"Constant analytical stress"
    )

    ax.text(
        0.02,
        0.06,
        stress_text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10.5,
        bbox=dict(
            boxstyle="round,pad=0.4",
            facecolor="#0d1117",
            edgecolor="#30363d",
        ),
    )

    plt.tight_layout()

    plt.savefig(
        results_dir / "stress_distribution.png",
        dpi=300,
        bbox_inches="tight",
        facecolor=fig.get_facecolor(),
    )

    plt.close(fig)

    print("Professional plots saved to:", results_dir)


if __name__ == "__main__":
    main()
