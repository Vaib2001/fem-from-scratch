from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from tapered_bar import (
    analytical_displacement,
    analytical_stress,
    solve_tapered_bar,
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
    })


def main():

    apply_github_dark_theme()

    # ---------------------------------------------------------
    # Engineering problem
    # ---------------------------------------------------------

    length = 1.0
    area_start = 0.01
    youngs_modulus = 210e9
    force = 100000.0

    # A(L) = 2 * A(0)
    taper_ratio = 1.0

    # Representative mesh for field plots
    num_elements_plot = 8

    # Meshes used for convergence study
    element_counts = np.array([
        1,
        2,
        4,
        8,
        16,
        32,
        64,
    ])

    # ---------------------------------------------------------
    # Analytical tip displacement
    # ---------------------------------------------------------

    exact_tip_displacement = analytical_displacement(
        length,
        length,
        area_start,
        youngs_modulus,
        force,
        taper_ratio,
    )

    # ---------------------------------------------------------
    # Mesh convergence study
    # ---------------------------------------------------------

    relative_errors = []

    fem_tip_displacements = []

    for num_elements in element_counts:

        (
            nodes,
            displacements,
            element_centres,
            stresses,
        ) = solve_tapered_bar(
            length,
            area_start,
            youngs_modulus,
            force,
            taper_ratio,
            int(num_elements),
        )

        fem_tip = displacements[-1]

        error = abs(
            fem_tip
            - exact_tip_displacement
        ) / abs(
            exact_tip_displacement
        )

        fem_tip_displacements.append(
            fem_tip
        )

        relative_errors.append(
            error
        )

    relative_errors = np.array(
        relative_errors
    )

    fem_tip_displacements = np.array(
        fem_tip_displacements
    )

    # ---------------------------------------------------------
    # Estimate convergence order
    #
    # error ~ N^(-p)
    # ---------------------------------------------------------

    slope, intercept = np.polyfit(
        np.log(element_counts),
        np.log(relative_errors),
        1,
    )

    convergence_order = -slope

    # ---------------------------------------------------------
    # Representative tapered-bar solution
    # ---------------------------------------------------------

    (
        nodes,
        displacements,
        element_centres,
        stresses,
    ) = solve_tapered_bar(
        length,
        area_start,
        youngs_modulus,
        force,
        taper_ratio,
        num_elements_plot,
    )

    x_exact = np.linspace(
        0.0,
        length,
        400,
    )

    u_exact = analytical_displacement(
        x_exact,
        length,
        area_start,
        youngs_modulus,
        force,
        taper_ratio,
    )

    sigma_exact = analytical_stress(
        x_exact,
        length,
        area_start,
        force,
        taper_ratio,
    )

    # ---------------------------------------------------------
    # Results directory
    # ---------------------------------------------------------

    project_root = (
        Path(__file__)
        .resolve()
        .parents[1]
    )

    results_dir = (
        project_root
        / "results"
    )

    results_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ---------------------------------------------------------
    # Plot 1:
    # Tapered-bar displacement
    # ---------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(9, 5.5)
    )

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
        color="#f0883e",
        marker="o",
        markersize=7,
        markeredgecolor="#f0f6fc",
        markeredgewidth=0.8,
        linewidth=1.8,
        label=f"FEM ({num_elements_plot} elements)",
    )

    ax.set_title(
        "Tapered Bar: FEM vs Analytical Displacement",
        pad=14,
    )

    ax.set_xlabel(
        "Position along bar [m]"
    )

    ax.set_ylabel(
        "Axial displacement [mm]"
    )

    ax.grid(True)

    ax.legend(
        loc="upper left"
    )

    plt.tight_layout()

    plt.savefig(
        results_dir
        / "tapered_displacement.png",
        dpi=300,
        bbox_inches="tight",
        facecolor=fig.get_facecolor(),
    )

    plt.close(fig)

    # ---------------------------------------------------------
    # Plot 2:
    # Tapered-bar stress
    # ---------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(9, 5.5)
    )

    ax.plot(
        x_exact,
        sigma_exact / 1e6,
        color="#58a6ff",
        linewidth=2.8,
        label="Analytical stress",
    )

    ax.plot(
        element_centres,
        stresses / 1e6,
        linestyle="None",
        marker="o",
        markersize=7,
        markerfacecolor="#3fb950",
        markeredgecolor="#f0f6fc",
        markeredgewidth=0.8,
        label="FEM element stress",
    )

    ax.set_title(
        "Tapered Bar: Axial Stress Distribution",
        pad=14,
    )

    ax.set_xlabel(
        "Position along bar [m]"
    )

    ax.set_ylabel(
        "Axial stress [MPa]"
    )

    ax.grid(True)

    ax.legend()

    plt.tight_layout()

    plt.savefig(
        results_dir
        / "tapered_stress.png",
        dpi=300,
        bbox_inches="tight",
        facecolor=fig.get_facecolor(),
    )

    plt.close(fig)

    # ---------------------------------------------------------
    # Plot 3:
    # Mesh convergence
    # ---------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(9, 5.5)
    )

    ax.loglog(
        element_counts,
        relative_errors,
        color="#f0883e",
        linewidth=2.4,
        marker="o",
        markersize=7,
        markeredgecolor="#f0f6fc",
        markeredgewidth=0.8,
        label="FEM relative error",
    )

    # Second-order reference slope
    reference_error = (
        relative_errors[0]
        * (
            element_counts[0]
            / element_counts
        ) ** 2
    )

    ax.loglog(
        element_counts,
        reference_error,
        color="#58a6ff",
        linestyle="--",
        linewidth=2.0,
        label=r"Reference: $O(N^{-2})$",
    )

    ax.set_title(
        "Mesh Convergence: Tapered Axial Bar",
        pad=14,
    )

    ax.set_xlabel(
        "Number of finite elements"
    )

    ax.set_ylabel(
        "Relative error in tip displacement"
    )

    ax.grid(
        True,
        which="both",
    )

    ax.legend()

    convergence_text = (
        f"Estimated order ≈ "
        f"{convergence_order:.2f}"
    )

    ax.text(
        0.97,
        0.08,
        convergence_text,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        bbox=dict(
            boxstyle="round,pad=0.4",
            facecolor="#0d1117",
            edgecolor="#30363d",
        ),
    )

    plt.tight_layout()

    plt.savefig(
        results_dir
        / "tapered_convergence.png",
        dpi=300,
        bbox_inches="tight",
        facecolor=fig.get_facecolor(),
    )

    plt.close(fig)

    # ---------------------------------------------------------
    # Print convergence table
    # ---------------------------------------------------------

    print("=" * 76)

    print(
        "TAPERED BAR - MESH CONVERGENCE STUDY"
    )

    print("=" * 76)

    print(
        f"\nAnalytical tip displacement: "
        f"{exact_tip_displacement:.8e} m"
    )

    print(
        f"Estimated convergence order: "
        f"{convergence_order:.3f}"
    )

    print(
        "\nElements   FEM tip displacement [m]"
        "   Relative error"
    )

    print("-" * 76)

    for (
        num_elements,
        fem_tip,
        error,
    ) in zip(
        element_counts,
        fem_tip_displacements,
        relative_errors,
    ):

        print(
            f"{num_elements:8d}"
            f"   {fem_tip:.8e}"
            f"              {error:.8e}"
        )

    print(
        "\nPlots saved to:",
        results_dir,
    )


if __name__ == "__main__":
    main()
