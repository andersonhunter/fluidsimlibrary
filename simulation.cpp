#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "fluid_formulae.h"

#ifndef SIZE
#define SIZE 16
#endif

#ifndef DEBUG
#define DEBUG true
#endif

// Eventually, add in ability to switch between bilinear interp and bicubic for performance vs accuracy
bool hifi = true;  // Determine if the simulation is higher fidelity

struct Point {
    // A point within the grid which holds different quantities
    float temperature;  // Temperature
    float vx;           // X-velocity
    float vy;           // Y-Velocity
    float pressure;     // Temporary pressure, should be 0 after projection
    float density;      // Density
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
    // Corrected Value = Forward Estimate + 0.5 * (Original Value - Backward Estimate)
    // Save the corrected point in the temporary struct, then swap out when done
    
    //TODO: Maybe add a minmod limiter or filtering to prevent oscillations
    // Currently doing bicubic interpolation for the 4x4 grid [i - 1, j - 1] to [i + 2, j + 2]
    // Bilinear would be more performant but less accurate, add as an option later

    // Create struct of temporary values
    struct Point* tempGrid = (struct Point *)malloc((SIZE * SIZE) * sizeof(Point));

    for(int x = 0; x < SIZE; x++) {
        // Need to add boundary checking
        //fprintf(stdout, "%d,,", x);
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
            // Maybe turn backward into a backStruct at some point?

            float xPrev = (float)x - forwardX * timestep;
            float yPrev = (float)y - forwardY * timestep;
            // Round to nearest whole index
            i = floor(xPrev);
            j = floor(yPrev);
            // Calculate fractional offsets
            dx = xPrev - i;
            dy = yPrev - j;

            float tempXPrev[4];
            float tempYPrev[4];

            for(int n = -1; n <= 2; n++) {
                tempXPrev[n + 1] = cubicInterpolate(
                    getSafeVx(i + n, j - 1, grid),
                    getSafeVx(i + n, j, grid),
                    getSafeVx(i + n, j + 1, grid),
                    getSafeVx(i + n, j + 2, grid),
                    dx
                );
                tempYPrev[n + 1] = cubicInterpolate(
                    getSafeVy(i + n, j - 1, grid),
                    getSafeVy(i + n, j, grid),
                    getSafeVy(i + n, j + 1, grid),
                    getSafeVy(i + n, j + 2, grid),
                    dx
                );
            }

            // Calculate backward estimates
            backwardX = cubicInterpolate(tempXPrev[0], tempXPrev[1], tempXPrev[2], tempXPrev[3], dy);
            backwardY = cubicInterpolate(tempYPrev[0], tempYPrev[1], tempYPrev[2], tempYPrev[3], dy);

            // Calculate corrected value
            //  Corrected Value = Forward Estimate + 0.5 * (Original Value - Backward Estimate)
            float correctedX = forwardX + 0.5 * (getAtIndex(x, y, grid).vx - backwardX);
            float correctedY = forwardY + 0.5 * (getAtIndex(x, y, grid).vy - backwardY);

            // Copy new value into temporary struct
            setVxAtIndex(x, y, correctedX, tempGrid);
            setVyAtIndex(x, y, correctedY, tempGrid);
            if(DEBUG)
            {    // Debug stuff
                //float oldX = getAtIndex(x, y, grid).vx;
                //float oldY = getAtIndex(x, y, grid).vy;
                //float oldMag = sqrt(oldX * oldX + oldY * oldY);
                //float newMag = sqrt(correctedX * correctedX + correctedY * correctedY);
                //fprintf(stdout, "%.3f,", newMag);
                }
        }
        if(DEBUG) {
            //fprintf(stdout, "\n")'
        }
    }

    // Project temporary grid onto permanent grid
    memcpy(grid, tempGrid, SIZE * SIZE * sizeof(Point));
    free(tempGrid);
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
    // Returns 0 if the index is OOB or within boundary layer
    return (x > 0 && x < (SIZE - 1) && y > 0 && y < (SIZE - 1)) ? getAtIndex(x, y, grid).vx : 0.;
}

float getSafeVy(int x, int y, struct Point* grid) {
    // Helper function to safely get the y-velocity at a given index
    // Returns 0 if the index is OOB or within boundary layer
    return (x > 0 && x < (SIZE - 1) && y > 0 && y < (SIZE - 1)) ? getAtIndex(x, y, grid).vy : 0.;
}

int main(int argc, char* argv[]) {
    
    // Initialize grid of points
    // All boundary points should have 0 velocity for no-slip
    // Pseudo-2D array, index in with pointer arithmetic grid[SIZE * row + column]
    // Currently setting up with a linear spread
    struct Point* grid = (struct Point *)malloc((SIZE * SIZE) * sizeof(Point));
    for (int row = 0; row < SIZE; row++) {
        for (int col = 0; col < SIZE; col++) {
            if (row == 0 || row == SIZE - 1 || col % (SIZE - 1) == 0) {
                grid[SIZE * row + col].temperature = 0.;
                grid[SIZE * row + col].vx          = 0.;
                grid[SIZE * row + col].vy          = 0.;
                grid[SIZE * row + col].pressure    = 0.;
                grid[SIZE * row + col].density     = 1.;
            }
            else {
                grid[SIZE * row + col].temperature = 0.;
                grid[SIZE * row + col].vx          = 1. - ((float)col / 15.);
                grid[SIZE * row + col].vy          = 1. - ((float)col / 15.);
                grid[SIZE * row + col].pressure    = 0.;
                grid[SIZE * row + col].density     = 1.;
            }
        }
    }
    if (DEBUG) {
        //fprintf(stdout, "0,1,2,3,4,5,6,7,8,9,10\n\n");
        fprintf(stdout, "t,x,y\n0,1,1,,2,2\n");
    }
    
    float xpos = 1.;
    float tempx = 1.;
    float ypos = 1.;
    float tempy = 1.;

    float xpos1 = 2.;
    float tempx1 = 2.;
    float ypos1 = 2.;
    float tempy1 = 2.;

    for(float i = 0.; i <= 15.1; i += 0.1) {
        //fprintf(stdout, "\ntimestep = %.1f\n", i);
        calculateAdvection(grid, i);

        // Debug stuff, track pathlines for two points
        if (DEBUG) {
            xpos += getAtIndex(floor(tempx), floor(tempy), grid).vx;
            ypos += getAtIndex(floor(tempx), floor(tempy), grid).vy;
            tempx = xpos;
            tempy = ypos;

            xpos1 += getAtIndex(floor(tempx1), floor(tempy1), grid).vx;
            ypos1 += getAtIndex(floor(tempx1), floor(tempy1), grid).vy;
            tempx1 = xpos1;
            tempy1 = ypos1;

            fprintf(stdout, "%.1f,%.3f,%.3f,,%.3f,%.3f\n", i, xpos, ypos, xpos1, ypos1);
        }
        
    }

    free(grid);
    return 0;
}
