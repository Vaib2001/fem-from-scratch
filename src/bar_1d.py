import numpy as np


def validate_inputs(length, area, youngs_modulus, force, num_elements):
    """
    Validate the physical and numerical input parameters.
    """

    if length <= 0:
        raise ValueError("Bar length must be greater than zero.")

    if area <= 0:
        raise ValueError("Cross-sectional area must be greater than zero.")

    if youngs_modulus <= 0:
        raise ValueError("Young's modulus must be greater than zero.")

    if not np.isfinite(force):
        raise ValueError("Applied force must be finite.")

    if not isinstance(num_elements, int) or num_elements < 1:
        raise ValueError("Number of elements must be a positive integer.")


def element_stiffness(youngs_modulus, area, element_length):
    """
    Return the stiffness matrix of a 2-node linear bar element.

    Parameters
    ----------
    youngs_modulus : float
        Young's modulus [Pa].
    area : float
        Cross-sectional area [m^2].
    element_length : float
        Length of the finite element [m].

    Returns
    -------
    ndarray
        2 x 2 element stiffness matrix.
    """

    return (
        youngs_modulus * area / element_length
    ) * np.array([
        [1.0, -1.0],
        [-1.0, 1.0]
    ])


def analytical_displacement(
    x,
    area,
    youngs_modulus,
    force
):
    """
    Analytical displacement field for a uniform axial bar.

    u(x) = P*x / (E*A)
    """

    x = np.asarray(x, dtype=float)

    return (
        force * x
        / (youngs_modulus * area)
    )


def analytical_stress(area, force):
    """
    Analytical axial stress for a uniform bar.

    sigma = P / A
    """

    return force / area


def solve_bar(
    length,
    area,
    youngs_modulus,
    force,
    num_elements
):
    """
    Solve a fixed-free 1D axial bar using the Finite Element Method.

    The bar is fixed at x = 0 and subjected to an axial point load
    at x = L.

    Parameters
    ----------
    length : float
        Total bar length [m].
    area : float
        Cross-sectional area [m^2].
    youngs_modulus : float
        Young's modulus [Pa].
    force : float
        Axial force applied at the free end [N].
    num_elements : int
        Number of finite elements.

    Returns
    -------
    nodes : ndarray
        Nodal coordinates [m].
    displacements : ndarray
        Nodal displacements [m].
    stresses : ndarray
        Element stresses [Pa].
    """

    validate_inputs(
        length,
        area,
        youngs_modulus,
        force,
        num_elements
    )

    # ---------------------------------------------------------
    # 1. Generate the finite-element mesh
    # ---------------------------------------------------------

    num_nodes = num_elements + 1

    nodes = np.linspace(
        0.0,
        length,
        num_nodes
    )

    element_length = length / num_elements

    # ---------------------------------------------------------
    # 2. Initialize global system
    #
    # K u = F
    # ---------------------------------------------------------

    global_stiffness = np.zeros(
        (num_nodes, num_nodes)
    )

    global_force = np.zeros(num_nodes)

    # ---------------------------------------------------------
    # 3. Compute element stiffness matrix
    # ---------------------------------------------------------

    ke = element_stiffness(
        youngs_modulus,
        area,
        element_length
    )

    # ---------------------------------------------------------
    # 4. Assemble global stiffness matrix
    # ---------------------------------------------------------

    for element in range(num_elements):

        element_dofs = np.array([
            element,
            element + 1
        ])

        global_stiffness[
            np.ix_(
                element_dofs,
                element_dofs
            )
        ] += ke

    # ---------------------------------------------------------
    # 5. Apply external point load
    # ---------------------------------------------------------

    global_force[-1] = force

    # ---------------------------------------------------------
    # 6. Apply essential boundary condition
    #
    # u(0) = 0
    # ---------------------------------------------------------

    free_dofs = np.arange(
        1,
        num_nodes
    )

    reduced_stiffness = global_stiffness[
        np.ix_(
            free_dofs,
            free_dofs
        )
    ]

    reduced_force = global_force[
        free_dofs
    ]

    # ---------------------------------------------------------
    # 7. Solve reduced finite-element system
    # ---------------------------------------------------------

    displacements = np.zeros(num_nodes)

    displacements[
        free_dofs
    ] = np.linalg.solve(
        reduced_stiffness,
        reduced_force
    )

    # ---------------------------------------------------------
    # 8. Recover element strains and stresses
    # ---------------------------------------------------------

    strains = (
        np.diff(displacements)
        / element_length
    )

    stresses = (
        youngs_modulus
        * strains
    )

    return nodes, displacements, stresses


if __name__ == "__main__":

    # ---------------------------------------------------------
    # Example engineering problem
    # ---------------------------------------------------------

    LENGTH = 1.0               # m
    AREA = 0.01                # m^2
    YOUNGS_MODULUS = 210e9     # Pa
    FORCE = 100000.0           # N
    NUM_ELEMENTS = 4

    nodes, displacements, stresses = solve_bar(
        LENGTH,
        AREA,
        YOUNGS_MODULUS,
        FORCE,
        NUM_ELEMENTS
    )

    # ---------------------------------------------------------
    # Analytical solution
    # ---------------------------------------------------------

    exact_displacements = analytical_displacement(
        nodes,
        AREA,
        YOUNGS_MODULUS,
        FORCE
    )

    exact_tip_displacement = analytical_displacement(
        LENGTH,
        AREA,
        YOUNGS_MODULUS,
        FORCE
    )

    exact_stress = analytical_stress(
        AREA,
        FORCE
    )

    absolute_error = abs(
        displacements[-1]
        - exact_tip_displacement
    )

    if abs(exact_tip_displacement) > 0:

        relative_error = (
            absolute_error
            / abs(exact_tip_displacement)
        )

    else:

        relative_error = 0.0

    # ---------------------------------------------------------
    # Display results
    # ---------------------------------------------------------

    print("=" * 55)
    print("1D AXIAL BAR - FINITE ELEMENT ANALYSIS")
    print("=" * 55)

    print(f"\nNumber of elements : {NUM_ELEMENTS}")
    print(f"Number of nodes    : {len(nodes)}")

    print("\nNode coordinates [m]")
    print(nodes)

    print("\nNodal displacements [m]")
    print(displacements)

    print("\nElement stresses [Pa]")
    print(stresses)

    print("\n" + "-" * 55)
    print("ANALYTICAL VALIDATION")
    print("-" * 55)

    print(
        f"FEM tip displacement        : "
        f"{displacements[-1]:.8e} m"
    )

    print(
        f"Analytical tip displacement : "
        f"{exact_tip_displacement:.8e} m"
    )

    print(
        f"Absolute error              : "
        f"{absolute_error:.8e} m"
    )

    print(
        f"Relative error              : "
        f"{relative_error:.8e}"
    )

    print(
        f"FEM stress                  : "
        f"{stresses[-1] / 1e6:.6f} MPa"
    )

    print(
        f"Analytical stress           : "
        f"{exact_stress / 1e6:.6f} MPa"
    )
