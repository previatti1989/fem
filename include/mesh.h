#ifndef MESH_H
#define MESH_H

typedef struct {
    int num_nodes;
    int num_elements;
    double* nodes;       // (x, y) coordinates stored as [x0, y0, x1, y1, ...]
    int* elements;       // Node indices for each element
    int* boundary_nodes; // Indices of boundary nodes (Dirichlet)
    int num_boundary;
} FEMMesh;

// Function prototypes
void generate_mesh(FEMMesh* mesh, int Nx, int Ny);
void free_mesh(FEMMesh* mesh);

#endif // MESH_H