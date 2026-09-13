import numpy as np

from src.bar_1d import solve_bar
from src.tapered_bar import (
    analytical_displacement,
    solve_tapered_bar,
)


def test_uniform_bar_fixed_boundary():
    """The fixed end must have zero displacement."""

    _, displacements, _ = solve_bar(
        length=1.0,
        area=0.01,
        youngs_modulus=210e9,
        force=100000.0,
        num_elements=4,
    )

    assert displacements[0] == 0.0


def test_uniform_bar_tip_displacement():
    """FEM tip displacement must match the analytical solution."""

    length = 1.0
    area = 0.01
    youngs_modulus = 210e9
    force = 100000.0

    _, displacements, _ = solve_bar(
        length,
        area,
        youngs_modulus,
        force,
        4,
    )

    exact = (
        force
        * length
        / (youngs_modulus * area)
    )

    assert np.isclose(
        displacements[-1],
        exact,
        rtol=1e-12,
    )


def test_uniform_bar_stress():
    """Uniform bar stress must equal P/A."""

    area = 0.01
    force = 100000.0

    _, _, stresses = solve_bar(
        length=1.0,
        area=area,
        youngs_modulus=210e9,
        force=force,
        num_elements=4,
    )

    exact_stress = force / area

    assert np.allclose(
        stresses,
        exact_stress,
        rtol=1e-12,
    )


def tapered_tip_error(num_elements):
    """Return relative error in tapered-bar tip displacement."""

    length = 1.0
    area_start = 0.01
    youngs_modulus = 210e9
    force = 100000.0
    taper_ratio = 1.0

    _, displacements, _, _ = solve_tapered_bar(
        length,
        area_start,
        youngs_modulus,
        force,
        taper_ratio,
        num_elements,
    )

    exact = analytical_displacement(
        length,
        length,
        area_start,
        youngs_modulus,
        force,
        taper_ratio,
    )

    return abs(
        displacements[-1] - exact
    ) / abs(exact)


def test_tapered_bar_converges():
    """Mesh refinement must reduce tapered-bar discretization error."""

    error_2 = tapered_tip_error(2)
    error_8 = tapered_tip_error(8)
    error_32 = tapered_tip_error(32)

    assert error_8 < error_2
    assert error_32 < error_8


def test_tapered_bar_fine_mesh_accuracy():
    """A sufficiently refined mesh should achieve high accuracy."""

    error = tapered_tip_error(32)

    assert error < 1e-4
