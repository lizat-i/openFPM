#pragma once
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
// A constant to indicate boundary particles
#define BOUNDARY 0
// A constant to indicate fluid particles
#define FLUID 1
// initial spacing between particles dp in the formulas
// const double dp = 0.04;
const double dp = 0.0085;
// Maximum height of the fluid water
// is going to be calculated and filled later on
double h_swl = 0.0;
// c_s in the formulas (constant used to calculate the sound speed)
const double coeff_sound = 20.0;
// gamma in the formulas
const double gamma_ = 7.0;

// // sqrt(3.0*dp*dp) support of the kernel
// const double H = 0.0147224318643;

// sqrt(3.0*dp*dp) support of the kernel
const double H = std::sqrt(3.0*dp*dp)   ;  
// Reference densitu 1000Kg/m^3
const double rho_zero = 1000.0;

// Eta in the formulas
const double Eta2 = 0.01 * H*H;
// alpha in the formula
const double visco = 0.1;


// // Mass of the fluid particles
// const double MassFluid = 0.000614125;
// // Mass of the boundary particles
// const double MassBound = 0.000614125;

// Mass of the fluid particles
const double initVolume = pow(dp,3);
// Mass of the boundary particles
const double MassBound = rho_zero*initVolume;
// Mass of the fluid particles
const double MassFluid = rho_zero*initVolume;
 
// End simulation time
#ifdef TEST_RUN
const double t_end = 0.001;
#else
const double t_end = 1.5;
#endif
// Gravity acceleration
const double gravity = 9.81;

// Filled later require h_swl, it is b in the formulas
double B = 0.0;
// Constant used to define time integration
const double timesteppingConstant = 0.2;
// Minimum T
const double DtMin = 0.00001;
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
// calculated force
const int force = 5;
// velocity
const int velocity = 6;
// velocity at previous step
const int velocity_prev = 7;
// velocity at previous step
const int x_pre = 8;

// velocity at previous step
const int viscousFoces = 9;

// velocity at previous step
const double bodyforce[3] = {gravity, 0.0, 0.0};


// Type of the vector containing particles
typedef vector_dist<3,double,aggregate<size_t,double,  double,    double,     double,     double[3], double[3], double[3]  , double[3]      , double[3]      >> particles;
//                                       |      |        |          |            |            |         |            |          |                |                
//                                       |      |        |          |            |            |         |            |          |                |         
//                                     type   density   density    Pressure    delta       force     velocity  velocity_prev position_prev    viscous_forces
//                                             

// NON Vanilla Entries

const double kinematic_viscosity = 1e-1;

// 3d Cavity 
double Xmin = 0.0;
double Ymin = 0.0 - (dp * 3.0);
double Zmin = 0.0;
double Xmax = 1.0;
double Ymax = 1.0 + (dp * 3.0);
double Zmax = 1.0;

double L_x = Xmax - Xmin;
double L_y = Ymax - Ymin;
double L_z = Zmax - Zmin;

double characteristic_Length = 1;

size_t NumberOfParticles_x = ceil(L_x / dp);
size_t NumberOfParticles_y = ceil(L_y / dp);
size_t NumberOfParticles_z = ceil(L_z / dp);


const double theoMaxVel = ( gravity * characteristic_Length*characteristic_Length )/ (2 *kinematic_viscosity/ rho_zero ) ; 


const double cbar_1 = sqrt((gravity * characteristic_Length) * (1/0.01)); 
const double cbar_2 = sqrt(theoMaxVel* theoMaxVel * (1/0.01)); 
const double cbar_3 = sqrt((theoMaxVel* kinematic_viscosity) / ( (1/0.01)*characteristic_Length)   ); 

const double cbar  = std::max({cbar_1, cbar_2, cbar_3})  ;


double eps = std::numeric_limits<double>::min();
// dynamic_viscosity in the formula

