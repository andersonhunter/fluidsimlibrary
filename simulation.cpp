#include <stdlib.h>
#include <stdio.h>
#include "fluid_formulae.h"

#ifndef SIZE
#define SIZE 16
#endif

struct Point {
    // A point within the grid which holds different quantities
    float pressure;
    float temperature;
};

void calculateAdvection(struct Point* grid) {

}

struct Point getAtIndex(int x, int y, struct Point* grid) {
    // Retrieve a point at the given index from the grid
    // Grid[x][y] = grid[SIZE * x + y]
    return grid[SIZE * x + y];
}

int main(int argc, char* argv[]) {
    
    // Initialize grid of points
    // Pseudo-2D array, index in with pointer arithmetic grid[SIZE * row + column]
    // Ex: grid[2][15] = grid[SIZE * 2 + 15] = grid[47]
    struct Point* grid = (struct Point *)malloc((SIZE * SIZE) * sizeof(Point));

    // Set initial values
    for(int i = 0; i < SIZE * SIZE; i++) {
        grid[i].temperature = 0.;
        grid[i].pressure = 0.;
    }

    free(grid);
    return 0;
}
