#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mesh.h"
#include "matrix_ops.h"
#include "elliptic_solver.h"
#include "solvers.h"

double k_function(double x, double y) {
    return 1.0;
}

int nonlinear_scalar = 1.0;

//  Linear Solution: u_exact = sin(pi x) sin(pi y)
double u_exact_dirichlet(double x, double y) {
    return sin(M_PI * x) * sin(M_PI * y);
}

//  Source term for Dirichlet BC: f = 2pi² sin(pix) sin(piy)
double f_function_dirichlet(double x, double y) {
    return 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
}

// Dirichlet BC:
double g_D_dirichlet(double x, double y) {
    return u_exact_dirichlet(x, y);
}

//  Linear Solution: u_exact = x(1-x)y(1-y)
double u_exact_mixed(double x, double y) {
    return x * (1 - x) * y * (1 - y);
}

//  Source term for Mixed BC: f = -2[y(1-y) + x(1-x)]
double f_function_mixed(double x, double y) {
    return -2 * (y * (1 - y) + x * (1 - x));
}

//  Neumann BC
double g_N_mixed(double x, double y) {
    return (y == 0) ? -2 * x * (1 - x) : (y == 1) ? 2 * x * (1 - x) : NAN;
}

//  Dirichlet BC
double g_D_mixed(double x, double y) {
    return (x == 0 || x == 1) ? 0.0 : NAN;
}

//  Nonlinear Poisson Solution: u_exact = sin(pix) sin(piy)
double u_exact_nonlinear(double x, double y) {
    return sin(M_PI * x) * sin(M_PI * y);
}

//  Source term for Nonlinear Poisson: f = 2pi² sin(pix) sin(piy) - k sin³(pix) sin³(piy)
double f_nonlinear_function(double x, double y) {
    double u = sin(M_PI * x) * sin(M_PI * y);
    return 2 * M_PI * M_PI * u - 1.0 * pow(u, 3);  // Assume scalar = 1.0
}

// Nonlinear Dirichlet BC
double g_D_nonlinear(double x, double y) {
    return 0.0;
}

double g_D_mixed_nonlinear(double x, double y) {
    return (x == 0 || x == 1) ? 0.0 : NAN; // Only apply Dirichlet on x-boundaries
}

// Nonlinear Neumann BC
double g_N_mixed_nonlinear(double x, double y) {
    if (x == 0 || x == 1) return NAN;  // Dirichlet boundaries
    if (y == 0) return M_PI * sin(M_PI * x);  // Bottom boundary
    if (y == 1) return -M_PI * sin(M_PI * x); // Top boundary
    return NAN;
}

// Compare numerical solution with exact analytical solution
void validate_numerical_solution(const FEMMesh* mesh, const FEMVector* numerical_solution, double (*u_exact)(double, double)) {
    printf("\nComparing FEM solution with analytical solution:\n");

    double error_norm = 0.0;
    double max_error = 0.0;

    for (int i = 0; i < mesh->num_nodes; i++) {
        double x = mesh->nodes[2 * i];
        double y = mesh->nodes[2 * i + 1];
        double exact_value = u_exact(x, y);
        double computed_value = numerical_solution->values[i];  // Placeholder
        double error = fabs(computed_value - exact_value);

        error_norm += error * error;
        if (error > max_error) max_error = error;

        printf("Node %d: u_exact = %.6f, u_FEM = %.6f, error = %.6f\n", i, exact_value, computed_value, error);
    }

    error_norm = sqrt(error_norm / mesh->num_nodes);
    printf("\nL2 Norm of Error: %.6e\n", error_norm);
    printf("Max Error: %.6e\n", max_error);
}

void solve_and_validate(FEMMatrix* A, FEMVector* b, FEMVector* u, const FEMMesh* mesh, const char* method, double (*u_exact)(double, double)) {
    initialize_vector(u, mesh->num_nodes);
    set_vector_zero(u);

    if (strcmp(method, "QR") == 0) {
        printf("\nSolving using QR Decomposition...\n");
        qr_solver(A, b, u);
    }
    else if (strcmp(method, "CG") == 0) {
        printf("\nSolving using Conjugate Gradient (CG)...\n");
        cg_solver(A, b, u, 1e-6, 1000);
    }
    else if (strcmp(method, "LSQR") == 0) {
        printf("\nSolving using LSQR...\n");
        lsqr_solver(A, b, u, 1e-6, 1000);
    }
    else if (strcmp(method, "GMRES") == 0) {
        printf("\nSolving using GMRES...\n");
        int k_max = A->rows;
        printf("\nSolving GMRES k_max = %d\n", k_max);
        gmres_solver(A, b, u, 1e-6, 1000, k_max);
    }

    validate_numerical_solution(mesh, u, u_exact);
    free_vector(u);
}

int main() {
    int Nx = 4, Ny = 4;
    const char* solver_method = "CG";

    FEMMesh mesh;
    EllipticFEMSystem fem;
    FEMVector u_solution;

    // Generate mesh
    generate_mesh(&mesh, Nx, Ny);

    //  Dirichlet Poisson Test
    printf("\n==============================\n");
    printf("Testing LINEAR Poisson Equation (Dirichlet BC):\n");
    printf("==============================\n");
    initialize_elliptic_system(&fem, mesh.num_nodes, k_function, f_function_dirichlet, LINEAR, DIRICHLET_ONLY, 0.0);
    assemble_elliptic_stiffness(&fem, &mesh);
    assemble_elliptic_load_vector(&fem, &mesh);
    apply_elliptic_boundary_conditions(&fem, &mesh, g_D_dirichlet, NULL);

    solve_and_validate(&fem.stiffness_matrix, &fem.load_vector, &u_solution, &mesh, solver_method, u_exact_dirichlet);
    free_elliptic_system(&fem);

    //  Mixed Poisson Test
    printf("\n==============================\n");
    printf("Testing LINEAR Poisson Equation (Mixed BC):\n");
    printf("==============================\n");
    initialize_elliptic_system(&fem, mesh.num_nodes, k_function, f_function_mixed, LINEAR, MIXED_BOUNDARY, 0.0);
    assemble_elliptic_stiffness(&fem, &mesh);
    assemble_elliptic_load_vector(&fem, &mesh);
    apply_elliptic_boundary_conditions(&fem, &mesh, g_D_mixed, g_N_mixed);

    solve_and_validate(&fem.stiffness_matrix, &fem.load_vector, &u_solution, &mesh, solver_method, u_exact_mixed);
    free_elliptic_system(&fem);

    //  Nonlinear Poisson Test
    printf("\n==============================\n");
    printf("Testing NONLINEAR Poisson Equation:\n");
    printf("==============================\n");
    initialize_elliptic_system(&fem, mesh.num_nodes, k_function, f_nonlinear_function, NONLINEAR, DIRICHLET_ONLY, nonlinear_scalar);
    assemble_elliptic_stiffness(&fem, &mesh);
    assemble_mass_matrix(&fem, &mesh);
    assemble_elliptic_load_vector(&fem, &mesh);
    apply_elliptic_boundary_conditions(&fem, &mesh, g_D_nonlinear, NULL);

    solve_and_validate(&fem.stiffness_matrix, &fem.load_vector, &u_solution, &mesh, solver_method, u_exact_nonlinear);
    free_elliptic_system(&fem);

    //  Nonlinear Poisson Test
    printf("\n==============================\n");
    printf("Testing NONLINEAR Poisson Equation:\n");
    printf("==============================\n");
    initialize_elliptic_system(&fem, mesh.num_nodes, k_function, f_nonlinear_function, NONLINEAR, MIXED_BOUNDARY, nonlinear_scalar);
    assemble_elliptic_stiffness(&fem, &mesh);
    assemble_mass_matrix(&fem, &mesh);
    assemble_elliptic_load_vector(&fem, &mesh);
    apply_elliptic_boundary_conditions(&fem, &mesh, g_D_mixed_nonlinear, g_N_mixed_nonlinear);

    solve_and_validate(&fem.stiffness_matrix, &fem.load_vector, &u_solution, &mesh, solver_method, u_exact_nonlinear);
    free_elliptic_system(&fem);

}