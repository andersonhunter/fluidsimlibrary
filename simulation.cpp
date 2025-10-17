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

int main(int argc, char* argv[]) {
    
    // Initialize grid of points
    struct Point* grid = (struct Point *)malloc(SIZE * sizeof(Point));

    // Set initial values
    for(int i = 0; i < SIZE; i++) {
        grid[i].pressure = 0;
        grid[i].temperature = 0;
    }

    return 0;
}
