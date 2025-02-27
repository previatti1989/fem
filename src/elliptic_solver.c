#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "elliptic_solver.h"
#include "matrix_ops.h"

// Initialize the system with equation type, boundary type, and nonlinear scalar
void initialize_elliptic_system(EllipticFEMSystem* fem, int num_nodes, 
    CoeffFunction k_func, SourceFunction f_func, EquationType eq_type, BoundaryType b_type, double nonlinear_scalar) {
    
    fem->num_nodes = num_nodes;
    initialize_matrix(&fem->stiffness_matrix, num_nodes, num_nodes);
    initialize_matrix(&fem->mass_matrix, num_nodes, num_nodes);
    initialize_vector(&fem->load_vector, num_nodes);
    initialize_vector(&fem->solution, num_nodes);

    fem->k_function = k_func;
    fem->f_function = f_func;
    fem->equation_type = eq_type;
    fem->boundary_type = b_type;
    fem->nonlinear_scalar = nonlinear_scalar;
}

// Assemble stiffness matrix for both LINEAR and NONLINEAR equations
void assemble_elliptic_stiffness(EllipticFEMSystem* fem, const FEMMesh* mesh) {
    for (int i = 0; i < mesh->num_elements; i++) {
        int n1 = mesh->elements[3 * i];
        int n2 = mesh->elements[3 * i + 1];
        int n3 = mesh->elements[3 * i + 2];

        double x1 = mesh->nodes[2 * n1], y1 = mesh->nodes[2 * n1 + 1];
        double x2 = mesh->nodes[2 * n2], y2 = mesh->nodes[2 * n2 + 1];
        double x3 = mesh->nodes[2 * n3], y3 = mesh->nodes[2 * n3 + 1];

        // Compute area of the triangle
        double area = 0.5 * fabs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));

        // Compute shape function derivatives
        double bvec[3] = { y2 - y3, y3 - y1, y1 - y2 };
        double cvec[3] = { x3 - x2, x1 - x3, x2 - x1 };

        // Compute average k value
        double k_avg = (fem->k_function(x1, y1) + fem->k_function(x2, y2) + fem->k_function(x3, y3)) / 3.0;

        // Evaluate k(x,y) at the element centroid
        double x_c = (x1 + x2 + x3) / 3.0;
        double y_c = (y1 + y2 + y3) / 3.0;
        double k_e = fem->k_function(x_c, y_c);  //  Correctly applying k(x,y)

        // Compute local stiffness matrix
        double Ke[3][3];
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                Ke[a][b] = k_e * (bvec[a] * bvec[b] + cvec[a] * cvec[b]) / (4.0 * area);
            }
        }

        // Assemble into global stiffness matrix
        int nodes[3] = { n1, n2, n3 };
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                fem->stiffness_matrix.values[nodes[a] * fem->num_nodes + nodes[b]] += Ke[a][b];
            }
        }
    }
}

// Assemble the mass matrix M
void assemble_mass_matrix(EllipticFEMSystem* fem, const FEMMesh* mesh) {
    for (int i = 0; i < mesh->num_elements; i++) {
        int n1 = mesh->elements[3 * i];
        int n2 = mesh->elements[3 * i + 1];
        int n3 = mesh->elements[3 * i + 2];

        double x1 = mesh->nodes[2 * n1], y1 = mesh->nodes[2 * n1 + 1];
        double x2 = mesh->nodes[2 * n2], y2 = mesh->nodes[2 * n2 + 1];
        double x3 = mesh->nodes[2 * n3], y3 = mesh->nodes[2 * n3 + 1];

        // Compute area of the triangle
        double area = 0.5 * fabs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));

        // Mass matrix for a triangular element (Lumped Mass Approximation)
        double Me[3][3] = {
            { area / 6.0, area / 12.0, area / 12.0 },
            { area / 12.0, area / 6.0, area / 12.0 },
            { area / 12.0, area / 12.0, area / 6.0 }
        };

        int nodes[3] = { n1, n2, n3 };
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                fem->mass_matrix.values[nodes[a] * fem->num_nodes + nodes[b]] += Me[a][b];
            }
        }
    }
}

// Apply nonlinear modification to stiffness matrix
void apply_nonlinear_stiffness(EllipticFEMSystem* fem, const FEMVector* u) {
    for (int i = 0; i < fem->num_nodes; i++) {
        double u3 = pow(u->values[i], 3);  // Compute u³
        for (int j = 0; j < fem->num_nodes; j++) {
            fem->stiffness_matrix.values[i * fem->num_nodes + j] +=
                fem->nonlinear_scalar * fem->mass_matrix.values[i * fem->num_nodes + j] * u3;
        }
    }
}



// Assemble load vector F
void assemble_elliptic_load_vector(EllipticFEMSystem* fem, const FEMMesh* mesh) {
#pragma omp parallel for
    for (int i = 0; i < mesh->num_nodes; i++) {
        double x = mesh->nodes[2 * i];
        double y = mesh->nodes[2 * i + 1];
        fem->load_vector.values[i] = fem->f_function(x, y);
    }
}

// Apply selected boundary conditions
void apply_elliptic_boundary_conditions(EllipticFEMSystem* fem, const FEMMesh* mesh, double (*g_D)(double, double), double (*g_N)(double, double)) {
#pragma omp parallel for
    for (int i = 0; i < mesh->num_boundary; i++) {
        int b_node = mesh->boundary_nodes[i];
        double x = mesh->nodes[2 * b_node];
        double y = mesh->nodes[2 * b_node + 1];

        if (fem->boundary_type == DIRICHLET_ONLY || fem->boundary_type == MIXED_BOUNDARY) {
            for (int j = 0; j < fem->num_nodes; j++) {
                fem->stiffness_matrix.values[b_node * fem->num_nodes + j] = 0.0;
                fem->stiffness_matrix.values[j * fem->num_nodes + b_node] = 0.0;
            }
            fem->stiffness_matrix.values[b_node * fem->num_nodes + b_node] = 1.0;
            fem->load_vector.values[b_node] = g_D(x,y);
        }

        if (fem->boundary_type == NEUMANN_ONLY || fem->boundary_type == MIXED_BOUNDARY) {
            double ds = 1.0 / (mesh->num_boundary);  // Approximate segment length
            fem->load_vector.values[b_node] += g_N(x, y) * ds;
        }
    }
}

// Free memory
void free_elliptic_system(EllipticFEMSystem* fem) {
    free_matrix(&fem->stiffness_matrix);
    free_matrix(&fem->mass_matrix);
    free_vector(&fem->load_vector);
    free_vector(&fem->solution);
}
