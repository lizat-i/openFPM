#pragma once
// #define SE_CLASS1
// #define STOP_ON_ERROR
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include <limits>
#include <iostream>
#include <fstream>

std::ofstream outFile;
// A constant to indicate boundary particles
#define BOUNDARY 0
// A constant to indicate fluid particles
#define FLUID 1
const double dp = 0.02;

const double EPS = std::numeric_limits<double>::epsilon() ; 
// Maximum height of the fluid water
// is going to be calculated and filled later on
double h_swl = 0.0;
// c_s in the formulas (constant used to calculate the sound speed)
const double coeff_sound = 20.0;
// gamma in the formulas
const double gamma_ = 7.0;
// sqrt(3.0*dp*dp) support of the kernel
const double H = sqrt(3.0 * dp * dp);
// number of kappa
const double kappa = 2.0;
// number of kappa
const double smoothingRadius = kappa * H;
// Eta in the formulas
const double Eta2 = 0.01 * H * H;
// alpha in the formula
const double visco = 0.1;
// cbar in the formula (calculated later)
double cbar = 0.0;
// Reference densitu 1000Kg/m^3
const double rho_zero = 1000.0;
// Mass of the fluid particles
const double MassFluid = std::pow(dp, 3) * rho_zero;
// Mass of the boundary particles
const double MassBound = std::pow(dp, 3) * rho_zero;
// End simulation time
#ifdef TEST_RUN
const double t_end = 0.001;
#else
const double t_end = 4.0;
#endif
// Gravity acceleration
const double gravity = 9.81;

// Filled later require h_swl, it is b in the formulas
double B = 0.0;
// Constant used to define time integration
const double delta_t_coeeficient = 0.2;
// Minimum T
const double DtMin = 1.00 * 1e-5;
// Minimum Rho allowed
const double RhoMin = 700.0;
// Maximum Rho allowed
const double RhoMax = 1300.0;
// Filled in initialization
double max_fluid_height = 0.0;
// Properties
// FLUID or BOUNDARY
const size_t type = 0;
// Density
const int rho = 1;
// Density at step n-1
const int rho_prev = 2;
// Pressure
const int Pressure = 3;
// Delta rho calculated in the force calculation
const int drho = 4;
// calculated pressure Force
const int pressure_acc = 5;
// velocity
const int velocity = 6;
// velocity at previous step
const int velocity_prev = 7;
// calculated viscous Force
const int viscous_acc = 8;
// position at previous timestep
const int x_pre = 9;
// position at previous timestep
const int rho_err = 10;
// position at previous timestep
const int Pressure_prev = 10;
// Gravity or external forces
const double bodyforce[3] = {0, 0, -gravity};

// maximum allowed density error
const double maxDensityVariation = 0.01;

// smoothing parameter
const double intPConst = 0.5;

double Xmin = 0.0;
double Ymin = 0.0;
double Zmin = 0.0;

double Xmax = 0.0;
double Ymax = 0.0;
double Zmax = 0.0;

double L_x = 0.0;
double L_y = 0.0;
double L_z = 0.0;