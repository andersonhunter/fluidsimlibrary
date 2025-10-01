#ifndef FLUID_FORMULAE_H_
#define FLUID_FORMULAE_H_

// Define constants
extern const float GRAVITY;   // Acceleration due to gravity in m/sec^2
extern const float SW_WATER;  // Specific weight of water at 4C
extern const float F_PI;      // Float PI

float forceMA(float mass, float acceleration);

float weightMA(float mass);

float absolutePressure(float atmosphericPressure, float vacuumPressure);

float density(float mass, float volume);

float specificWeightWV(float weight, float volume);

float specificWeightDensity(float density);

float specificVolumeSpecificWeight(float specificWeight);

float specificVolumeDensity(float density);

float specificGravity(float density);

float surfaceTension(float attraction, float length);

float vaporSolidInterfaceSurfaceTension(
    float liquidSolidSurfaceTension,
    float vaporLiquidSurfaceTension,
    float contactAngle
);

float columnarWeight(float fluidSpecificWeight, float height, float diameter);

float forceSurfaceTension(float surfaceTension, float diameter);

float columnHeight(
    float surfaceTension,
    float contactAngle,
    float specificWeight,
    float diameter
);

float stress(float force, float area);

float shearStress(float viscosityCoefficient, float velocity, float distance);

float stokesLaw(float outsideRadius, float viscosity, float sphereSW, float fluidSW);

float kinematicViscosity(float viscosity, float density);

float pressureFA(float force, float area);

#endif