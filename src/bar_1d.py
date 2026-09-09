import numpy as np


def element_stiffness(E, A, length):
    """
    Return the stiffness matrix of a 2-node 1D bar element.
    """
    return (E * A / length) * np.array([
        [1.0, -1.0],
        [-1.0, 1.0]
    ])


def solve_bar(length, area, youngs_modulus, force, num_elements):
    """
    Solve a fixed-free 1D axial bar using the Finite Element Method.

    Parameters
    ----------
    length : float
        Total bar length [m]
    area : float
        Cross-sectional area [m^2]
    youngs_modulus : float
        Young's modulus [Pa]
    force : float
        Axial force applied at the free end [N]
    num_elements : int
        Number of finite elements

    Returns
    -------
    nodes : ndarray
        Node coordinates
    displacements : ndarray
        Nodal displacements
    stresses : ndarray
        Stress in each element
    """

    num_nodes = num_elements + 1

    # ---------------------------------------------------------
    # 1. Generate mesh
    # ---------------------------------------------------------
    nodes = np.linspace(0.0, length, num_nodes)

    element_length = length / num_elements

    # ---------------------------------------------------------
    # 2. Initialize global stiffness matrix and force vector
    # ---------------------------------------------------------
    K = np.zeros((num_nodes, num_nodes))
    F = np.zeros(num_nodes)

    # ---------------------------------------------------------
    # 3. Assemble global stiffness matrix
    # ---------------------------------------------------------
    ke = element_stiffness(
        youngs_modulus,
        area,
        element_length
    )

    for element in range(num_elements):

        node_1 = element
        node_2 = element + 1

        element_nodes = [node_1, node_2]

        for i in range(2):
            for j in range(2):
                K[element_nodes[i], element_nodes[j]] += ke[i, j]

    # ---------------------------------------------------------
    # 4. Apply external force
    # ---------------------------------------------------------
    F[-1] = force

    # ---------------------------------------------------------
    # 5. Apply boundary condition
    #
    # Node 0 is fixed:
    # u(0) = 0
    # ---------------------------------------------------------
    free_dofs = np.arange(1, num_nodes)

    K_reduced = K[np.ix_(free_dofs, free_dofs)]
    F_reduced = F[free_dofs]

    # ---------------------------------------------------------
    # 6. Solve K u = F
    # ---------------------------------------------------------
    displacements = np.zeros(num_nodes)

    displacements[free_dofs] = np.linalg.solve(
        K_reduced,
        F_reduced
    )

    # ---------------------------------------------------------
    # 7. Calculate element stresses
    # ---------------------------------------------------------
    stresses = np.zeros(num_elements)

    for element in range(num_elements):

        u1 = displacements[element]
        u2 = displacements[element + 1]

        strain = (u2 - u1) / element_length

        stresses[element] = youngs_modulus * strain

    return nodes, displacements, stresses


if __name__ == "__main__":

    # Example engineering problem

    L = 1.0          # m
    A = 0.01         # m^2
    E = 210e9        # Pa
    P = 100000.0     # N

    number_of_elements = 4

    nodes, u, stress = solve_bar(
        L,
        A,
        E,
        P,
        number_of_elements
    )

    # Analytical displacement at free end
    analytical_displacement = P * L / (E * A)

    print("Node coordinates:")
    print(nodes)

    print("\nNodal displacements [m]:")
    print(u)

    print("\nElement stresses [Pa]:")
    print(stress)

    print("\nFEM displacement at free end:")
    print(u[-1])

    print("\nAnalytical displacement at free end:")
    print(analytical_displacement)

    relative_error = abs(
        u[-1] - analytical_displacement
    ) / analytical_displacement

    print("\nRelative error:")
    print(relative_error)
