#pragma once
#pragma once
// problem parameters
#ifndef CONFIG_HH
#define CONFIG_HH
using real = float;

const int     NX = 32;                  // number of node points along X and Y (cavity length in lattice units)
const int     NY = 32;
const int     NZ = 128;

const double  REYNOLDS_NUMBER = 1E6;    // REYNOLDS_NUMBER = LID_VELOCITY * N / kinematicViscosity

// don't change these unless you know what you are doing

const int     Qcc = 19;                    // number of discrete velocity aections used
const double  DENSITY = 2.7;            // fluid density in lattice units
const double  LID_VELOCITY = 0.07;      // lid velocity in lattice units

#endif