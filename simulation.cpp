#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fluid_formulae.h"

#ifndef SIZE
#define SIZE 16
#endif

bool hifi = true;  // Determine if the simulation is higher fidelity

struct Point {
    // A point within the grid which holds different quantities
    float temperature = 0.;  // Temperature
    float vx = 0.;           // X-velocity
    float vy = 0.;           // Y-Velocity
    float pressure = 0.;     // Temporary pressure, should be 0 after projection
    float density = 1.;      // Density
};

void calculateAdvection(struct Point* grid, float timestep) {
    // Use the MacCormack Method
    // Do a forward pass to estimate advection at each point using the Semi-Lagrangian method
    //   Trace backward along the velocity field to gain the origin point of the fluid
    //   Interpolate at that point
    //   Store as forward estimate
    // Do a backward pass using Semi-Lagrangian
    //   Trace forward to estimate where the fluid ends up
    //   Interpolate value at that point
    //   Store as backward estimate
    //  Corrected Value = Forward Estimate + 0.5 * (Original Value - Backward Estimate)
    
    //TODO: Maybe add a minmod limiter or filtering to prevent oscillations

    for(int x = 0; x < SIZE; x++) {
        // Need to add boundary checking
        for(int y = 0; y < SIZE; y++) {
            // Calculate the forward estimate
            // Get previous (x, y) coords
            float xPrev = (float)x - getAtIndex(x, y, grid).vx * timestep;
            float yPrev = (float)y - getAtIndex(x, y, grid).vy * timestep;
            // Calculate fractional offsets (amount each point differs from its nearest point)
            float dx = xPrev - floor(xPrev);
            float dy = yPrev - floor(yPrev);
            // Now need to interpolate the x and y velocities using the bilinear interpolation thang
        }
    }
}

struct Point getAtIndex(int x, int y, struct Point* grid) {
    // Retrieve a point at the given index from the grid
    // Grid[x][y] = grid[SIZE * x + y]
    return grid[SIZE * x + y];
}

int main(int argc, char* argv[]) {
    
    // Initialize grid of points
    // Pseudo-2D array, index in with pointer arithmetic grid[SIZE * row + column]
    struct Point* grid = (struct Point *)malloc((SIZE * SIZE) * sizeof(Point));

    free(grid);
    return 0;
}
