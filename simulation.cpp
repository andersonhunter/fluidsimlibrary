#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fluid_formulae.h"

#ifndef SIZE
#define SIZE 16
#endif

// Eventually, add in ability to switch between bilinear interp and bicubic for performance vs accuracy
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
float cubicInterpolate(float p0, float p1, float p2, float p3, float dt);
float getSafeVx(int x, int y, struct Point* grid);
float getSafeVy(int x, int y, struct Point* grid);

float cubicInterpolate(float p0, float p1, float p2, float p3, float dt) {
    // Perform Centripetal Catmull-Rom interpolation to determine advected value
    // Calculate spline values, use coefficients from Catmull-Rom spline matrix
    float a0 = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
    float a1 = p0 - 2.5 * p1 + 2. * p2 - 0.5 * p3;
    float a2 = -0.5 * p0 + 0.5 * p2;
    // Return the cubic polynomial
    return dt * (dt * (a0 * dt + a1) + a2) + p1;
}

void calculateAdvection(struct Point* grid, float timestep) {
    // Use the MacCormack Method
    // Do a forward pass to estimate advection at each point using the Semi-Lagrangian method
    //   Trace forward along the velocity field to the estimated new location of the fluid
    //   Interpolate at that point
    //   Store as forward estimate
    // Do a backward pass using Semi-Lagrangian
    //   Trace backward from the new point
    //   Interpolate value at that point
    //   Store as backward estimate
    //  Corrected Value = Forward Estimate + 0.5 * (Original Value - Backward Estimate)
    
    //TODO: Maybe add a minmod limiter or filtering to prevent oscillations
    // Currently doing bicubic interpolation for the 4x4 grid [i - 1, j - 1] to [i + 2, j + 2]
    // Bilinear would be more performant but less accurate, add as an option later
    for(int x = 0; x < SIZE; x++) {
        // Need to add boundary checking
        for(int y = 0; y < SIZE; y++) {
            float forwardX, forwardY;
            float backwardX, backwardY;
            // Calculate the forward estimate
            // Get next (x, y) coords
            float xNext = (float)x + getAtIndex(x, y, grid).vx * timestep;
            float yNext = (float)y + getAtIndex(x, y, grid).vy * timestep;
            // Round to nearest whole index
            int i = floor(xNext);
            int j = floor(yNext);
            // Calculate fractional offsets
            float dx = xNext - i;
            float dy = yNext - j;
            // Interpolate x and y velocities at the new point
            // vx, vy = CI(CI(p0..p3, dx), dy)
            // If point is outside the grid, use v = 0 for no-slip boundary conditions
            // Add more quantities (temperature, density, etc) here
            // Maybe turn forward into a forwardStruct at some point?

            // Estimate row-by-row x-velocity
            float tempX1, tempX2, tempX3, tempX4, tempY1, tempY2, tempY3, tempY4;
            float tempX[4];
            float tempY[4];
            for(int n = -1; n <= 2; n++) {
                tempX[n + 1] = cubicInterpolate(
                    getSafeVx(i + n, j - 1, grid),
                    getSafeVx(i + n, j, grid),
                    getSafeVx(i + n, j + 1, grid),
                    getSafeVx(i + n, j + 2, grid),
                    dx
                );
                tempY[n + 1] = cubicInterpolate(
                    getSafeVy(i + n, j - 1, grid),
                    getSafeVy(i + n, j, grid),
                    getSafeVy(i + n, j + 1, grid),
                    getSafeVy(i + n, j + 2, grid),
                    dx
                );
            }

            // Calculate forward estimates
            forwardX = cubicInterpolate(tempX[0], tempX[1], tempX[2], tempX[3], dy);
            forwardY = cubicInterpolate(tempY[0], tempY[1], tempY[2], tempY[3], dy);

            // Back trace the point and interpolate x and y velocities
            // vx, vy = CI(CI(p0..p3, dx), dy)
            // If point is outside the grid, use v = 0 for no-slip boundary conditions
            // Add more quantities (temperature, density, etc) here
            // Maybe turn forward into a forwardStruct at some point?

            float xPrev = (float)x - forwardX * timestep;
            float yPrev = (float)y - forwardY * timestep;
            // Round to nearest whole index
            i = floor(xPrev);
            j = floor(yPrev);
            // Calculate fractional offsets
            dx = xPrev - i;
            dy = yPrev - j;

            tempX1 = cubicInterpolate(
                getSafeVx(i - 1, j - 1, grid),
                getSafeVx(i - 1, j, grid),
                getSafeVx(i - 1, j + 1, grid),
                getSafeVx(i - 1, j + 2, grid),
                dx
            );
            tempX2 = cubicInterpolate(
                getSafeVx(i, j - 1, grid),
                getSafeVx(i, j, grid),
                getSafeVx(i, j + 1, grid),
                getSafeVx(i, j + 2, grid),
                dx
            );
            tempX3 = cubicInterpolate(
                getSafeVx(i + 1, j - 1, grid),
                getSafeVx(i + 1, j, grid),
                getSafeVx(i + 1, j + 1, grid),
                getSafeVx(i + 1, j + 2, grid),
                dx
            );
            tempX4 = cubicInterpolate(
                getSafeVx(i + 2, j - 1, grid),
                getSafeVx(i + 2, j, grid),
                getSafeVx(i + 2, j + 1, grid),
                getSafeVx(i + 2, j + 2, grid),
                dx
            );

            // Estimate row-by-row y-velocity
            tempY1 = cubicInterpolate(
                getSafeVy(i - 1, j - 1, grid),
                getSafeVy(i - 1, j, grid),
                getSafeVy(i - 1, j + 1, grid),
                getSafeVy(i - 1, j + 2, grid),
                dx
            );
            tempY2 = cubicInterpolate(
                getSafeVy(i, j - 1, grid),
                getSafeVy(i, j, grid),
                getSafeVy(i, j + 1, grid),
                getSafeVy(i, j + 2, grid),
                dx
            );
            tempY3 = cubicInterpolate(
                getSafeVy(i + 1, j - 1, grid),
                getSafeVy(i + 1, j, grid),
                getSafeVy(i + 1, j + 1, grid),
                getSafeVy(i + 1, j + 2, grid),
                dx
            );
            tempY4 = cubicInterpolate(
                getSafeVy(i + 2, j - 1, grid),
                getSafeVy(i + 2, j, grid),
                getSafeVy(i + 2, j + 1, grid),
                getSafeVy(i + 2, j + 2, grid),
                dx
            );

            // Calculate backward estimates
            backwardX = cubicInterpolate(tempX1, tempX2, tempX3, tempX4, dy);
            backwardY = cubicInterpolate(tempY1, tempY2, tempY3, tempY4, dy);

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

float getSafeVx(int x, int y, struct Point* grid) {
    // Helper function to safely get the x-velocity at a given index
    // Returns 0 if the index is OOB
    return (x >= 0 && x < SIZE && y >= 0 && y < SIZE) ? getAtIndex(x, y, grid).vx : 0.;
}

float getSafeVy(int x, int y, struct Point* grid) {
    // Helper function to safely get the y-velocity at a given index
    // Returns 0 if the index is OOB
    return (x >= 0 && x < SIZE && y >= 0 && y < SIZE) ? getAtIndex(x, y, grid).vy : 0.;
}

int main(int argc, char* argv[]) {
    
    // Initialize grid of points
    // Pseudo-2D array, index in with pointer arithmetic grid[SIZE * row + column]
    struct Point* grid = (struct Point *)malloc((SIZE * SIZE) * sizeof(Point));
    calculateAdvection(grid, 1);
    free(grid);
    return 0;
}
