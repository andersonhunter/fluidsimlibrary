#include <stdlib.h>
#include <stdio.h>
#include "fluid_formulae.h"

struct Point {
    float pressure = 0.;
    float velocity = 0.;
};

int main() {
    struct Point p1;
    printf("%f, %f", p1.pressure, p1.velocity);
    return 0;
}
