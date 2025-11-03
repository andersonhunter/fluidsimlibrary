#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
// Uncomment if doing any conversions or calculations
//#include "fluid_formulae.h"

#ifndef SIZE
#define SIZE 16
#endif

#ifndef CELLSIZE
#define CELLSIZE SIZE / SIZE
#endif

#ifndef DEBUG
#define DEBUG true
#endif

// Number of Jacobi iterations to perform
// More iterations = more resolution (but more computation time ):
#ifndef JACOBIS
#define JACOBIS 30
#endif

// Eventually, add in ability to switch between bilinear interp and bicubic for performance vs accuracy
bool hifi = true;  // Determine if the simulation is higher fidelity

struct Point {
    // A point within the grid which holds different quantities
    float temperature;  // Temperature
    float vx;           // X-velocity
    float vy;           // Y-Velocity
    float oldPressure;  // Previous pressure for Poisson
    float pressure;     // Temporary pressure, should be ~0 after projection
    float density;      // Density
    float sourceTerm;   // Divergence at the given point
};

struct Point* grid;

// Function prototypes
void calculateAdvection(float timestep);
struct Point getAtIndex(int x, int y);
void setVxAtIndex(int x, int y, float vx);
void setVyAtIndex(int x, int y, float vy);
float cubicInterpolate(float p0, float p1, float p2, float p3, float dt);
float getSafeVx(int x, int y);
float getSafeVy(int x, int y);
void initGrid();
float solveDivergence(float p0, float p1, float p2, float p3, float sourceTerm);
float getSafePressure(float* pressureGrid, int x, int y);

void initGrid() {
    // Initialize 2D grid with 0 values for all except density
    grid = (struct Point *)malloc((SIZE * SIZE) * sizeof(Point));
    for (int row = 0; row < SIZE; row++) {
        for (int col = 0; col < SIZE; col++) {
            grid[SIZE * row + col].temperature = 0.;
                grid[SIZE * row + col].vx          = 0.;
                grid[SIZE * row + col].vy          = 0.;
                grid[SIZE * row + col].pressure    = 0.;
                grid[SIZE * row + col].density     = 1.;
        }
    }
}

float solveDivergence(float p0, float p1, float p2, float p3, float sourceTerm) {
    // Solve Poisson for divergence at each point 
    return 0.25 * (p0 + p1 + p2 + p3 - sourceTerm);
}

float cubicInterpolate(float p0, float p1, float p2, float p3, float dt) {
    // Perform Centripetal Catmull-Rom interpolation to determine advected value
    // Calculate spline values, use coefficients from Catmull-Rom spline matrix
    float a0 = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
    float a1 = p0 - 2.5 * p1 + 2. * p2 - 0.5 * p3;
    float a2 = -0.5 * p0 + 0.5 * p2;
    // Return the cubic polynomial
    return dt * (dt * (a0 * dt + a1) + a2) + p1;
}

void calculateAdvection(float timestep) {
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
            float xNext = (float)x + getAtIndex(x, y).vx * timestep;
            float yNext = (float)y + getAtIndex(x, y).vy * timestep;
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
                    getSafeVx(i + n, j - 1),
                    getSafeVx(i + n, j),
                    getSafeVx(i + n, j + 1),
                    getSafeVx(i + n, j + 2),
                    dx
                );
                tempY[n + 1] = cubicInterpolate(
                    getSafeVy(i + n, j - 1),
                    getSafeVy(i + n, j),
                    getSafeVy(i + n, j + 1),
                    getSafeVy(i + n, j + 2),
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
                    getSafeVx(i + n, j - 1),
                    getSafeVx(i + n, j),
                    getSafeVx(i + n, j + 1),
                    getSafeVx(i + n, j + 2),
                    dx
                );
                tempYPrev[n + 1] = cubicInterpolate(
                    getSafeVy(i + n, j - 1),
                    getSafeVy(i + n, j),
                    getSafeVy(i + n, j + 1),
                    getSafeVy(i + n, j + 2),
                    dx
                );
            }

            // Calculate backward estimates
            backwardX = cubicInterpolate(tempXPrev[0], tempXPrev[1], tempXPrev[2], tempXPrev[3], dy);
            backwardY = cubicInterpolate(tempYPrev[0], tempYPrev[1], tempYPrev[2], tempYPrev[3], dy);

            // Calculate corrected value
            //  Corrected Value = Forward Estimate + 0.5 * (Original Value - Backward Estimate)
            float correctedX = forwardX + 0.5 * (getAtIndex(x, y).vx - backwardX);
            float correctedY = forwardY + 0.5 * (getAtIndex(x, y).vy - backwardY);

            // Copy new value into temporary struct
            tempGrid[SIZE * x + y].vx = correctedX;
            tempGrid[SIZE * x + y].vy = correctedY;
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

    // Solve for divergence (ensure no areas in the fluid are compressed)
    // All points should have initial pressures = 0
    // Loop k times (figure out what an acceptable divergence looks like and modify later)
    // Initialize pressures to 0 in preparation for Poisson-ification
    float *tempPressures = (float *)calloc(SIZE * SIZE, sizeof(float));
    float *tempOldPressures = (float *)calloc(SIZE * SIZE, sizeof(float));

    // Calculate sourceTerm ((x-comp + y-comp) * material derivative * h^2) for each point
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            tempGrid[SIZE * i + j].sourceTerm = (getSafeVx(i + 1, j) - getSafeVx(i - 1, j)) / (2. * (float)CELLSIZE);
            tempGrid[SIZE * i + j].sourceTerm += (getSafeVy(i, j + 1) - getSafeVy(i, j - 1)) / (2. * (float)CELLSIZE);
            tempGrid[SIZE * i + j].sourceTerm *= (tempGrid[SIZE * i + j].density / timestep);
            tempGrid[SIZE * i + j].sourceTerm *= ((float)(CELLSIZE * CELLSIZE));
        }
    }
    // Jacobi
    for (int k = 0; k < JACOBIS; k++) {
        // Move over old pressures
        memcpy(tempOldPressures, tempPressures, SIZE * SIZE * sizeof(float));
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                tempPressures[SIZE * i + j] = solveDivergence(
                    getSafePressure(tempOldPressures, i + 1, j),
                    getSafePressure(tempOldPressures, i - 1, j),
                    getSafePressure(tempOldPressures, i, j + 1),
                    getSafePressure(tempOldPressures, i, j - 1),
                    tempGrid[SIZE * i + j].sourceTerm
                );
            }
        }
    }

    // Now correct x- and y-velocities to account for pressure buildup (or lack thereof)
    // vx = (temp vx) - ((dt / rho) * (Pright - Pleft) / 2 * CELLSIZE)
    // vy = (temp vy) - ((dt / rho) * (Pabove - Pbelow) / 2 * CELLSIZE)

    for (int x = 0; x < SIZE; x++) {
        for (int y = 0; y < SIZE; y++) {
            float scale = timestep / grid[SIZE * x + y].density;
            float h = 2. * (float)CELLSIZE;
            float xcorrection = scale * ((getSafePressure(tempPressures, x + 1, y) - getSafePressure(tempPressures, x - 1, y)) / h);
            float ycorrection = scale * ((getSafePressure(tempPressures, x, y - 1) - getSafePressure(tempPressures, x, y + 1)) / h);
            tempGrid[x * SIZE + y].vx -= xcorrection;
            tempGrid[x * SIZE + y].vy -= ycorrection;
        }
    }

    // Project temporary grid onto permanent grid
    memcpy(grid, tempGrid, SIZE * SIZE * sizeof(Point));
    free(tempGrid);
    free(tempPressures);
    free(tempOldPressures);
}

struct Point getAtIndex(int x, int y) {
    // Retrieve a point at the given index from the grid
    // Grid[x][y] = grid[SIZE * x + y]
    return grid[SIZE * x + y];
}

void setVxAtIndex(int x, int y, float vx) {
    // Set x-velocity at a given index
    grid[SIZE * x + y].vx = vx;
}

void setVyAtIndex(int x, int y, float vy) {
    grid[SIZE * x + y].vy = vy;
}

float getSafeVx(int x, int y) {
    // Helper function to safely get the x-velocity at a given index
    // Returns 0 if the index is OOB or within boundary layer
    return (x > 0 && x < (SIZE - 1) && y > 0 && y < (SIZE - 1)) ? getAtIndex(x, y).vx : 0.;
}

float getSafeVy(int x, int y) {
    // Helper function to safely get the y-velocity at a given index
    // Returns 0 if the index is OOB or within boundary layer
    return (x > 0 && x < (SIZE - 1) && y > 0 && y < (SIZE - 1)) ? getAtIndex(x, y).vy : 0.;
}

float getSafePressure(float* pressureGrid, int x, int y) {
    return (x > 0 && x < (SIZE - 1) && y > 0 && y < (SIZE - 1)) ? pressureGrid[SIZE * x + y] : 0.;
}

int main(int argc, char* argv[]) {
    initGrid();
    for (int i = 0; i < 100; i++) {
        calculateAdvection((float)i);
    }
    free(grid);
    return 0;
}
