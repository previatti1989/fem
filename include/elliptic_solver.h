#ifndef ELLIPTIC_SOLVER_H
#define ELLIPTIC_SOLVER_H

#include "matrix_ops.h"
#include "mesh.h"

// Equation type: Linear or Nonlinear
typedef enum {
    LINEAR,
    NONLINEAR
} EquationType;

// Boundary condition type
typedef enum {
    DIRICHLET_ONLY,
    NEUMANN_ONLY,
    MIXED_BOUNDARY
} BoundaryType;

// Define function pointer types
typedef double (*CoeffFunction)(double, double);
typedef double (*SourceFunction)(double, double);

typedef struct {
    int num_nodes;
    FEMMatrix stiffness_matrix;
    FEMMatrix mass_matrix;
    FEMVector load_vector;
    FEMVector solution;
    CoeffFunction k_function;
    SourceFunction f_function;
    EquationType equation_type;
    BoundaryType boundary_type;
    double nonlinear_scalar; // Scalar coefficient for u^3 term (only used in NONLINEAR mode)
} EllipticFEMSystem;

// Initialize system
void initialize_elliptic_system(EllipticFEMSystem* fem, int num_nodes, CoeffFunction k_func, SourceFunction f_func, EquationType eq_type, BoundaryType b_type, double nonlinear_scalar);

// Assemble system
void assemble_elliptic_stiffness(EllipticFEMSystem* fem, const FEMMesh* mesh);
void assemble_mass_matrix(EllipticFEMSystem* fem, const FEMMesh* mesh);
void assemble_elliptic_load_vector(EllipticFEMSystem* fem, const FEMMesh* mesh);

// Apply selected boundary conditions
void apply_elliptic_boundary_conditions(EllipticFEMSystem* fem, const FEMMesh* mesh, double(*g_D)(double, double), double (*g_N)(double, double));

// Apply nonlinear modification
void apply_nonlinear_stiffness(EllipticFEMSystem* fem, const FEMVector* u);

// Free memory
void free_elliptic_system(EllipticFEMSystem* fem);

#endif // ELLIPTIC_SOLVER_H
