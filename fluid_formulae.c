#include "fluid_formulae.h"
#include <math.h>

// Define constants
extern const float GRAVITY = 9.806;   // Acceleration due to gravity in m/sec^2
extern const float SW_WATER;  // Specific weight of water at 4 Celsius
extern const float F_PI = 3.141592654;

float forceMA(float mass, float acceleration) {
    // Force an object exerts
    // Units = N = kg * m/sec^2

    // Sanity check
    if (acceleration < 0.) {
        acceleration = 0.;
    }

    return mass * acceleration;
}

float weightMA(float mass) {
    // Weight given mass
    // Units = N = kg * m/sec^2

    return mass * GRAVITY;
}

float pressureFA(float force, float area) {
    // Pressure given force and area
    // Units = Pa = N/m^2

    return force * area;
}

float absolutePressure(float atmosphericPressure, float vacuumPressure) {
    // Absolute pressure
    // Units = Pa = N/m^2

    return atmosphericPressure - vacuumPressure;
}

float density(float mass, float volume) {
    // Density
    // Units = kg/m^3

    return mass / volume;
}

float specificWeightWV(float weight, float volume) {
    // Specific weight given weight and volume
    // Units = N/m^3 = kg/(m^2 * s^2)

    return weight / volume;

}

float specificWeightDensity(float density) {
    // Specific weight given density
    // Units = N/m^3 = kg/(m^2 * s^2)

    return density * GRAVITY;
}

float specificVolumeSpecificWeight(float specificWeight) {
    // Specific volume given specific weight
    // Units = m^3/kg

    return 1 / specificWeight;

}

float specificVolumeDensity(float density) {
    // Specific volume given density
    // Units = m^3/kg

    return 1 / (density * GRAVITY);
}

float specificGravity(float density) {
    // Specific gravity given density
    // Units = Unitless

    return density / SW_WATER;
}

float surfaceTension(float attraction, float length) {
    // Surface tension of a fluid
    // Molecular force of attraction per unit length of free surface
    // Units = N/m

    return attraction / length;
}

float vaporSolidInterfaceSurfaceTension(
    float liquidSolidSurfaceTension,
    float vaporLiquidSurfaceTension,
    float contactAngle
) {
    // Surface tension at the vapor-solid interface
    // Units = N/m

    return liquidSolidSurfaceTension + vaporLiquidSurfaceTension * cosf(contactAngle);
}

float columnarWeight(float fluidSpecificWeight, float height, float diameter) {
    // Weight of a column (i.e. force exerted downwards) 
    // Units = N = kg * m/sec^2

    (fluidSpecificWeight * height * F_PI * diameter * diameter) / 4;
}

float forceSurfaceTension(float surfaceTension, float diameter) {
    // Total force on a column due to surface tension
    // Units = N = kg/(m^2 * s^2)

    return surfaceTension * F_PI * diameter;
}

float stress(float force, float area) {
    // Basic stress
    // Units = Pa = N/m^2

    return force / area;
}

float shearStress(float viscosityCoefficient, float velocity, float distance) {
    // Force of a substance slipping along another substance
    // Units = Pa = N/m^2

    return (viscosityCoefficient * velocity) / distance;

}

float stokesLaw(float outsideRadius, float viscosity, float sphereSW, float fluidSW) {
    // Terminal velocity of a sphere through a fluid
    // Units = m/s

    return (2 * outsideRadius * outsideRadius) / (viscosity * (sphereSW - fluidSW));
}

float kinematicViscosity(float viscosity, float density) {
    // Viscosity per unit of density
    // Units = stoke = m^2/sec

    return viscosity / density;
}