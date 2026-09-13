import numpy as np


def validate_inputs(
    length,
    area_start,
    youngs_modulus,
    force,
    taper_ratio,
    num_elements,
):
    """
    Validate physical and numerical parameters.
    """

    if length <= 0:
        raise ValueError("Bar length must be greater than zero.")

    if area_start <= 0:
        raise ValueError("Initial cross-sectional area must be positive.")

    if youngs_modulus <= 0:
        raise ValueError("Young's modulus must be positive.")

    if not np.isfinite(force):
        raise ValueError("Applied force must be finite.")

    if not isinstance(num_elements, int) or num_elements < 1:
        raise ValueError(
            "Number of elements must be a positive integer."
        )

    # A(L) = A0 * (1 + taper_ratio)
    if 1.0 + taper_ratio <= 0:
        raise ValueError(
            "Taper ratio produces a non-positive end area."
        )


def area_function(
    x,
    length,
    area_start,
    taper_ratio,
):
    """
    Linearly varying cross-sectional area.

    A(x) = A0 * (1 + alpha*x/L)

    where alpha is the taper ratio.
    """

    x = np.asarray(x, dtype=float)

    return area_start * (
        1.0
        + taper_ratio * x / length
    )


def element_stiffness(
    youngs_modulus,
    area_average,
    element_length,
):
    """
    Stiffness matrix for a 2-node axial bar element.

    For a linearly varying area, the exact element-integrated
    stiffness can be expressed using the average area.
    """

    return (
        youngs_modulus
        * area_average
        / element_length
    ) * np.array([
        [1.0, -1.0],
        [-1.0, 1.0],
    ])


def analytical_displacement(
    x,
    length,
    area_start,
    youngs_modulus,
    force,
    taper_ratio,
):
    """
    Analytical displacement field for the tapered bar.

    A(x) = A0 * (1 + alpha*x/L)

    For alpha != 0:

    u(x) =
        P*L/(E*A0*alpha)
        * ln(1 + alpha*x/L)
    """

    x = np.asarray(x, dtype=float)

    if abs(taper_ratio) < 1e-14:
        return (
            force
            * x
            / (
                youngs_modulus
                * area_start
            )
        )

    return (
        force
        * length
        / (
            youngs_modulus
            * area_start
            * taper_ratio
        )
        * np.log(
            1.0
            + taper_ratio
            * x
            / length
        )
    )


def analytical_stress(
    x,
    length,
    area_start,
    force,
    taper_ratio,
):
    """
    Exact axial stress distribution.

    sigma(x) = P / A(x)
    """

    area = area_function(
        x,
        length,
        area_start,
        taper_ratio,
    )

    return force / area


def solve_tapered_bar(
    length,
    area_start,
    youngs_modulus,
    force,
    taper_ratio,
    num_elements,
):
    """
    Solve a fixed-free tapered axial bar using linear
    finite elements.

    The cross-sectional area varies linearly:

        A(x) = A0 * (1 + alpha*x/L)

    Parameters
    ----------
    length : float
        Total bar length [m].

    area_start : float
        Cross-sectional area at x = 0 [m^2].

    youngs_modulus : float
        Young's modulus [Pa].

    force : float
        Axial end load [N].

    taper_ratio : float
        Relative change in area from x = 0 to x = L.

        Example:
        taper_ratio = 1.0

        gives

        A(L) = 2 * A0

    num_elements : int
        Number of finite elements.

    Returns
    -------
    nodes : ndarray
        Nodal coordinates [m].

    displacements : ndarray
        Nodal displacements [m].

    element_centres : ndarray
        Element-centre coordinates [m].

    stresses : ndarray
        FEM element stresses [Pa].
    """

    validate_inputs(
        length,
        area_start,
        youngs_modulus,
        force,
        taper_ratio,
        num_elements,
    )

    # ---------------------------------------------------------
    # 1. Generate mesh
    # ---------------------------------------------------------

    num_nodes = num_elements + 1

    nodes = np.linspace(
        0.0,
        length,
        num_nodes,
    )

    # ---------------------------------------------------------
    # 2. Initialize global FEM system
    # ---------------------------------------------------------

    global_stiffness = np.zeros(
        (num_nodes, num_nodes)
    )

    global_force = np.zeros(num_nodes)

    # ---------------------------------------------------------
    # 3. Element formulation and global assembly
    # ---------------------------------------------------------

    for element in range(num_elements):

        node_1 = element
        node_2 = element + 1

        x1 = nodes[node_1]
        x2 = nodes[node_2]

        element_length = x2 - x1

        area_1 = area_function(
            x1,
            length,
            area_start,
            taper_ratio,
        )

        area_2 = area_function(
            x2,
            length,
            area_start,
            taper_ratio,
        )

        # Exact average area for a linear A(x)
        area_average = 0.5 * (
            area_1 + area_2
        )

        ke = element_stiffness(
            youngs_modulus,
            area_average,
            element_length,
        )

        element_dofs = np.array([
            node_1,
            node_2,
        ])

        global_stiffness[
            np.ix_(
                element_dofs,
                element_dofs,
            )
        ] += ke

    # ---------------------------------------------------------
    # 4. Apply axial load at free end
    # ---------------------------------------------------------

    global_force[-1] = force

    # ---------------------------------------------------------
    # 5. Apply essential boundary condition
    #
    # u(0) = 0
    # ---------------------------------------------------------

    free_dofs = np.arange(
        1,
        num_nodes,
    )

    reduced_stiffness = global_stiffness[
        np.ix_(
            free_dofs,
            free_dofs,
        )
    ]

    reduced_force = global_force[
        free_dofs
    ]

    # ---------------------------------------------------------
    # 6. Solve reduced FEM system
    # ---------------------------------------------------------

    displacements = np.zeros(
        num_nodes
    )

    displacements[
        free_dofs
    ] = np.linalg.solve(
        reduced_stiffness,
        reduced_force,
    )

    # ---------------------------------------------------------
    # 7. Recover element stresses
    # ---------------------------------------------------------

    element_centres = 0.5 * (
        nodes[:-1]
        + nodes[1:]
    )

    element_lengths = np.diff(nodes)

    strains = (
        np.diff(displacements)
        / element_lengths
    )

    stresses = (
        youngs_modulus
        * strains
    )

    return (
        nodes,
        displacements,
        element_centres,
        stresses,
    )


if __name__ == "__main__":

    # ---------------------------------------------------------
    # Example tapered-bar problem
    # ---------------------------------------------------------

    LENGTH = 1.0
    AREA_START = 0.01
    YOUNGS_MODULUS = 210e9
    FORCE = 100000.0

    # End area becomes twice the initial area
    TAPER_RATIO = 1.0

    NUM_ELEMENTS = 8

    (
        nodes,
        displacements,
        element_centres,
        stresses,
    ) = solve_tapered_bar(
        LENGTH,
        AREA_START,
        YOUNGS_MODULUS,
        FORCE,
        TAPER_RATIO,
        NUM_ELEMENTS,
    )

    # ---------------------------------------------------------
    # Analytical validation
    # ---------------------------------------------------------

    exact_tip_displacement = (
        analytical_displacement(
            LENGTH,
            LENGTH,
            AREA_START,
            YOUNGS_MODULUS,
            FORCE,
            TAPER_RATIO,
        )
    )

    fem_tip_displacement = (
        displacements[-1]
    )

    relative_error = abs(
        fem_tip_displacement
        - exact_tip_displacement
    ) / abs(
        exact_tip_displacement
    )

    # ---------------------------------------------------------
    # Display results
    # ---------------------------------------------------------

    print("=" * 60)
    print("TAPERED 1D AXIAL BAR - FEM ANALYSIS")
    print("=" * 60)

    print(
        f"\nNumber of elements : "
        f"{NUM_ELEMENTS}"
    )

    print(
        f"Initial area       : "
        f"{AREA_START:.6f} m^2"
    )

    print(
        f"End area           : "
        f"{area_function(
            LENGTH,
            LENGTH,
            AREA_START,
            TAPER_RATIO
        ):.6f} m^2"
    )

    print(
        f"\nFEM tip displacement        : "
        f"{fem_tip_displacement:.8e} m"
    )

    print(
        f"Analytical tip displacement : "
        f"{exact_tip_displacement:.8e} m"
    )

    print(
        f"Relative error              : "
        f"{relative_error:.8e}"
    )

    print(
        "\nElement stresses [MPa]"
    )

    print(
        stresses / 1e6
    )
