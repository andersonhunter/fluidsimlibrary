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

// Function prototypes
void calculateAdvection(struct Point* grid, float timestep);
struct Point getAtIndex(int x, int y, struct Point* grid);
void setVxAtIndex(int x, int y, float vx, struct Point* grid);
void setVyAtIndex(int x, int y, float vy, struct Point* grid);

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
    // Currently doing bicubic interpolation for the 4x4 grid [i - 1, j - 1] to [i + 2, j + 2]
    // Bilinear would be more performant but less accurate
    for(int x = 0; x < SIZE; x++) {
        // Need to add boundary checking
        for(int y = 0; y < SIZE; y++) {
            float forwardX, forwardY;
            float backwardX, backwardY;
            // Calculate the forward estimate
            // Get previous (x, y) coords
            float xPrev = (float)x - getAtIndex(x, y, grid).vx * timestep;
            float yPrev = (float)y - getAtIndex(x, y, grid).vy * timestep;
            // Clamp to nearest valid index and sanity check
            int i = floor(xPrev);
            int j = floor(yPrev);
            if (i < 0) {
                i = 0;
            }
            if (i >= SIZE) {
                i = SIZE - 1;
            }
            if (j < 0) {
                j = 0;
            }
            if (j >= SIZE) {
                j = SIZE - 1;
            }
            // Calculate fractional offsets
            float dx = xPrev - i;
            float dy = yPrev - j;
            // Forward interpolate x and y velocities
            // Add more quantities (temperature, density, etc) here
            // Maybe turn forward into a forwardStruct at some point?
            
        }
    }
}

struct Point getAtIndex(int x, int y, struct Point* grid) {
    // Retrieve a point at the given index from the grid
    // Grid[x][y] = grid[SIZE * x + y]
    return grid[SIZE * x + y];
}

void setVxAtIndex(int x, int y, float vx, struct Point* grid) {
    // Set x-velocity at a given index
    grid[SIZE * x + y].vx = vx;
}

void setVyAtIndex(int x, int y, float vy, struct Point* grid) {
    grid[SIZE * x + y].vy = vy;
}

int main(int argc, char* argv[]) {
    
    // Initialize grid of points
    // Pseudo-2D array, index in with pointer arithmetic grid[SIZE * row + column]
    struct Point* grid = (struct Point *)malloc((SIZE * SIZE) * sizeof(Point));
    free(grid);
    return 0;
}
