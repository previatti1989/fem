#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

// Generate a structured Nx x Ny grid for the unit square
void generate_mesh(FEMMesh* mesh, int Nx, int Ny) {
    mesh->num_nodes = (Nx + 1) * (Ny + 1);
    mesh->num_elements = 2 * Nx * Ny; // Two triangles per quadrilateral
    mesh->nodes = (double*)malloc(2 * mesh->num_nodes * sizeof(double));
    mesh->elements = (int*)malloc(3 * mesh->num_elements * sizeof(int));
    mesh->boundary_nodes = (int*)malloc(4 * (Nx + Ny) * sizeof(int)); // Max possible boundary nodes

    int node_index = 0, element_index = 0, boundary_index = 0;

    // Generate nodes (structured grid)
#pragma omp parallel for collapse(2) reduction(+:boundary_index)
    for (int j = 0; j <= Ny; j++) {
        for (int i = 0; i <= Nx; i++) {
            double x = (double)i / Nx;
            double y = (double)j / Ny;
            mesh->nodes[2 * node_index] = x;
            mesh->nodes[2 * node_index + 1] = y;

            // Identify boundary nodes (left, right, top, bottom)
            if (i == 0 || i == Nx || j == 0 || j == Ny) {
                mesh->boundary_nodes[boundary_index++] = node_index;
            }

            node_index++;
        }
    }

    // Generate triangular elements (each quadrilateral = 2 triangles)
#pragma omp parallel for
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            int n1 = j * (Nx + 1) + i;
            int n2 = n1 + 1;
            int n3 = n1 + (Nx + 1);
            int n4 = n3 + 1;

            // Triangle 1: (n1, n2, n3)
            mesh->elements[element_index++] = n1;
            mesh->elements[element_index++] = n2;
            mesh->elements[element_index++] = n3;

            // Triangle 2: (n2, n4, n3)
            mesh->elements[element_index++] = n2;
            mesh->elements[element_index++] = n4;
            mesh->elements[element_index++] = n3;
        }
    }

    mesh->num_boundary = boundary_index;
}

// Free memory
void free_mesh(FEMMesh* mesh) {
    free(mesh->nodes);
    free(mesh->elements);
    free(mesh->boundary_nodes);
}