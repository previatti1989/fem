#include <stdio.h>
#include <math.h>
#include "mesh.h"
#include "matrix_ops.h"
#include "elliptic_solver.h"

// Analytical solution for verification
double u_exact(double x, double y) {
    return sin(M_PI * x) * sin(M_PI * y);
}

// Source function
double f_function(double x, double y) {
    return 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
}

// Nonlinear source function
double f_nonlinear_function(double x, double y) {
    double u = u_exact(x, y);
    double scalar = 1.0;
    return 2 * M_PI * M_PI * u - scalar * pow(u, 3);
}

// Constant coefficient k(x,y) = 1
double k_function(double x, double y) {
    return 1.0;
}

// Compare numerical solution with exact analytical solution
void validate_numerical_solution(const FEMMesh* mesh, const double* numerical_solution) {
    printf("\nComparing FEM solution with analytical solution:\n");

    double error_norm = 0.0;
    double max_error = 0.0;

    for (int i = 0; i < mesh->num_nodes; i++) {
        double x = mesh->nodes[2 * i];
        double y = mesh->nodes[2 * i + 1];
        double exact_value = u_exact(x, y);
        double computed_value = numerical_solution[i];  // Placeholder
        double error = fabs(computed_value - exact_value);

        error_norm += error * error;
        if (error > max_error) max_error = error;

        printf("Node %d: u_exact = %.6f, u_FEM = %.6f, error = %.6f\n", i, exact_value, computed_value, error);
    }

    error_norm = sqrt(error_norm / mesh->num_nodes);
    printf("\nL2 Norm of Error: %.6e\n", error_norm);
    printf("Max Error: %.6e\n", max_error);
}

int main() {
    int Nx = 4, Ny = 4;
    FEMMesh mesh;
    EllipticFEMSystem fem;

    // Generate mesh
    generate_mesh(&mesh, Nx, Ny);

    // Test for LINEAR Poisson equation
    printf("\n==============================\n");
    printf("Testing LINEAR Poisson Equation:\n");
    printf("==============================\n");
    initialize_elliptic_system(&fem, mesh.num_nodes, k_function, f_function, LINEAR, DIRICHLET_ONLY, 0.0);
    assemble_elliptic_stiffness(&fem, &mesh);
    assemble_elliptic_load_vector(&fem, &mesh);
    apply_elliptic_boundary_conditions(&fem, &mesh);

    // Print stiffness matrix and load vector
    print_matrix(&fem.stiffness_matrix, "K (Linear)");
    print_vector(&fem.load_vector, "F (Linear)");

    validate_numerical_solution(&mesh, fem.load_vector.values);  // Placeholder for numerical solution
    free_elliptic_system(&fem);

    // Test for NONLINEAR Poisson equation
    printf("\n==============================\n");
    printf("Testing NONLINEAR Poisson Equation:\n");
    printf("==============================\n");
    initialize_elliptic_system(&fem, mesh.num_nodes, k_function, f_nonlinear_function, NONLINEAR, DIRICHLET_ONLY, 1.0);
    assemble_elliptic_stiffness(&fem, &mesh);
    assemble_mass_matrix(&fem, &mesh);
    assemble_elliptic_load_vector(&fem, &mesh);
    apply_elliptic_boundary_conditions(&fem, &mesh);

    // Placeholder initial guess for U (can be improved)
    FEMVector u_guess;
    initialize_vector(&u_guess, mesh.num_nodes);
    set_vector_zero(&u_guess);

    // Apply nonlinear stiffness modification
    apply_nonlinear_stiffness(&fem, &u_guess);

    // Print matrices for debugging
    print_matrix(&fem.stiffness_matrix, "K (Nonlinear Modified)");
    print_matrix(&fem.mass_matrix, "M (Nonlinear)");
    print_vector(&fem.load_vector, "F (Nonlinear)");

    validate_numerical_solution(&mesh, fem.load_vector.values);  // Placeholder for numerical solution
    free_elliptic_system(&fem);

    free_vector(&u_guess);
    free_elliptic_system(&fem);
    free_mesh(&mesh);
    return 0;
}