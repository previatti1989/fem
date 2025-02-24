#include <stdio.h>
#include "mesh.h"  // No change needed if compiling inside /FEM/

// Print mesh information
void print_mesh(const FEMMesh* mesh) {
    printf("Nodes (%d total):\n", mesh->num_nodes);
    for (int i = 0; i < mesh->num_nodes; i++) {
        printf("Node %d: (%.2f, %.2f)\n", i, mesh->nodes[2 * i], mesh->nodes[2 * i + 1]);
    }

    printf("\nElements (%d total):\n", mesh->num_elements);
    for (int i = 0; i < mesh->num_elements; i++) {
        printf("Element %d: [%d, %d, %d]\n", i, mesh->elements[3 * i], mesh->elements[3 * i + 1], mesh->elements[3 * i + 2]);
    }

    printf("\nBoundary Nodes (%d total):\n");
    for (int i = 0; i < mesh->num_boundary; i++) {
        printf("%d ", mesh->boundary_nodes[i]);
    }
    printf("\n");
}

int main() {
    int Nx = 4, Ny = 4; // Grid size
    FEMMesh mesh;
    
    generate_mesh(&mesh, Nx, Ny);
    
    print_mesh(&mesh);

    free_mesh(&mesh);
    return 0;
}