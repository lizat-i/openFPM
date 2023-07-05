// #define SE_CLASS1
// #define STOP_ON_ERROR
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include <limits>
#include <iostream>
#include <fstream>
#include "constants.h"

template <typename CellList>
inline void Init_Loop(particles &vd, CellList &NN, double &max_visc)
{
    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);
    while (part.isNext())
    {
        // ... a
        auto a = part.get();

        // Reset the Pressure force
        vd.template getProp<pressure_acc>(a)[0] = 0.0;
        vd.template getProp<pressure_acc>(a)[1] = 0.0;
        vd.template getProp<pressure_acc>(a)[2] = 0.0;

        vd.template getProp<viscous_acc>(a)[0] = 0.0;
        vd.template getProp<viscous_acc>(a)[1] = 0.0;
        vd.template getProp<viscous_acc>(a)[2] = 0.0;

        vd.template getProp<Pressure>(a) = Pbackground;

        Point<3, double> xa = vd.getPos(a);
        Point<3, double> va = vd.getProp<velocity>(a);

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);

        // We threat FLUID particle differently from BOUNDARY PARTICLES ...
        if (vd.getProp<type>(a) == FLUID)
        {

            auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
            double viscousFocesTermx = 0.0;
            double viscousFocesTermy = 0.0;
            double viscousFocesTermz = 0.0;
            double kernel = 0;

            // For each neighborhood particle
            while (Np.isNext() == true)
            {
                // ... q
                auto b = Np.get();

                Point<3, double> xb = vd.getPos(b);
                Point<3, double> vb = vd.getProp<velocity>(b);
                Point<3, double> dr = xa - xb;

                double rhob = vd.getProp<rho>(b);
                double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
                double r2 = norm2(dr);
                double r = sqrt(r2);

                // If the particles interact ...
                if (r < smoothingRadius)
                {

                    Point<3, double> v_rel = va - vb;
                    Point<3, double> DW;
                    DWab(dr, DW, r, false);
                    double wab = Wab(r);
                    kernel += wab;

                    double viscousForcesFactor = -massb * Pi(dr, r2, v_rel, rhoa, rhob, massb, max_visc);
                    viscousFocesTermx += viscousForcesFactor * DW.get(0);
                    viscousFocesTermy += viscousForcesFactor * DW.get(1);
                    viscousFocesTermz += viscousForcesFactor * DW.get(2);
                }
                ++Np;
            }

            vd.getProp<viscous_acc>(a)[0] = (kernel > EPS) ? viscousFocesTermx : 0.0;
            vd.getProp<viscous_acc>(a)[1] = (kernel > EPS) ? viscousFocesTermy : 0.0;
            vd.getProp<viscous_acc>(a)[2] = (kernel > EPS) ? viscousFocesTermz : 0.0;
        }
        ++part;
    }
}

template <typename CellList>
void position_and_velocity_prediction(particles &vd, CellList &NN, double dt)
{
    auto part = vd.getDomainIterator();

    vd.map();
    vd.ghost_get<type, rho, Pressure, Pressure_prev, velocity, velocity_prev, pressure_acc, viscous_acc>();
    vd.updateCellList(NN);

    while (part.isNext())
    {
        // ... a
        auto a = part.get();
        // skip Particle if it is of type Boundary
        if (vd.template getProp<type>(a) == BOUNDARY)
        {
            ++part;
            continue;
        }

        // predict velocity

        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + (vd.template getProp<pressure_acc>(a)[0] + vd.template getProp<viscous_acc>(a)[0] + bodyforce[0]) * dt;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + (vd.template getProp<pressure_acc>(a)[1] + vd.template getProp<viscous_acc>(a)[1] + bodyforce[1]) * dt;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + (vd.template getProp<pressure_acc>(a)[2] + vd.template getProp<viscous_acc>(a)[2] + bodyforce[2]) * dt;

        // predict position

        vd.getPos(a)[0] = vd.template getProp<x_pre>(a)[0] + 0.5 * (vd.template getProp<velocity_prev>(a)[0] + vd.template getProp<velocity>(a)[0]) * dt;
        vd.getPos(a)[1] = vd.template getProp<x_pre>(a)[1] + 0.5 * (vd.template getProp<velocity_prev>(a)[1] + vd.template getProp<velocity>(a)[1]) * dt;
        vd.getPos(a)[2] = vd.template getProp<x_pre>(a)[2] + 0.5 * (vd.template getProp<velocity_prev>(a)[2] + vd.template getProp<velocity>(a)[2]) * dt;

        ++part;

        /*
        If particle leaves domain during prediction step
        remove it from parallelized vector to avoid memorycorruption
         and print message
        */

        if (vd.getPos(a)[0] < Xmin || vd.getPos(a)[1] < Ymin || vd.getPos(a)[2] < Zmin ||
            vd.getPos(a)[0] > Xmax || vd.getPos(a)[1] > Ymax || vd.getPos(a)[2] > Zmax)
        {
            to_remove.add(a.getKey());
            outFile << "added to remove :  " << std::endl;
        }
    }

    vd.remove(to_remove, 0);
    cnt++;
}

template <typename CellList>
inline void predictDensity(particles &vd, CellList &NN, double &max_visc, const double &dt, double &errorGlobal)
{
    /*
     paralelization housekeeping
     currently irrelevant
     because testing with one processor
   */
    vd.map();
    vd.ghost_get<type, rho, Pressure, Pressure_prev, velocity, velocity_prev, pressure_acc, viscous_acc>();
    vd.updateCellList(NN);

    auto part = vd.getDomainIterator();

    while (part.isNext())
    {
        auto a = part.get();

        // skip boundaryParticles
        if (vd.getProp<type>(a) == BOUNDARY)
        {
            ++part;
            continue;
        };

        // Get particle properties a
        Point<3, double> xa = vd.getPos(a);
        Point<3, double> va = vd.getProp<velocity>(a);

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);
        double Pa = vd.getProp<Pressure>(a);

        // Initialize quantities for summation in loop
        Point<3, double> gradientSum_1 = {0, 0, 0};
        Point<3, double> gradientSum = {0, 0, 0};
        Point<3, double> drScaled = {0, 0, 0};

        double kernel = 0;
        double densitySumation = 0.0;
        double densityDiffusionTerm = 0.0;
        double pressureNeighbouring = 0.0;
        double sumSquaredGradient = 0.0;

        vd.getProp<drho>(a) = 0.0;

        // get Neighborhood of particle
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            // ... q
            auto b = Np.get();

            // Get particle properties b

            Point<3, double> xb = vd.getPos(b);
            Point<3, double> vb = vd.getProp<velocity>(b);

            double Pb = vd.getProp<Pressure>(b);
            double rhob = vd.getProp<rho>(b);
            double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;

            // calculate quantities
            Point<3, double> dr = xa - xb;
            Point<3, double> v_rel = va - vb;
            Point<3, double> DW;

            double r2 = norm2(dr);
            double r = sqrt(r2);
            DWab(dr, DW, r, false);
            double wab = Wab(r);

            // If the particles interact ...
            if (r < smoothingRadius)
            {
                // Density pediction
                // Desnsity sumation
                double dot_vrel_DW = dotProduct(v_rel, DW);
                densitySumation += massb * dot_vrel_DW;

                // Desnsity diffusion
                // drScaled = 2 * massb * (vd.getProp<rho_prev>(a) - vd.getProp<rho_prev>(b)) / (vd.getProp<rho_prev>(b) * (r2 + Eta2)) * dr;
                drScaled = 2 * massb * (rhoa - rhob) / (rhob * (r2 + Eta2)) * dr;
                densityDiffusionTerm += dotProduct(drScaled, DW);
            }
            ++Np;
        }

        // Density pediction
        vd.getProp<drho>(a) = densitySumation + visco * H * cbar * densityDiffusionTerm;

        // pressureLogger << "density shifting :         " << densitySumation << "\n";
        // pressureLogger << "density diffusion :    " << visco * H * cbar * densityDiffusionTerm << "\n";

        vd.getProp<rho>(a) = vd.getProp<rho_prev>(a) + vd.getProp<drho>(a) * dt;

        // Density error calculation and tracking

        double rho_error = vd.getProp<rho>(a) - rho_zero;
        double rho_error_normalized = rho_error / rho_zero;
        errorGlobal = std::max(errorGlobal, rho_error_normalized);
        vd.getProp<rho_err>(a) = rho_error;
        ++part;
    }
}

template <typename CellList>
inline void predictBeta(particles &vd, CellList &NN, const double &dt, double &maxShephard, double &beta)
{
    /*
     paralelization housekeeping
     currently irrelevant
     because testing with one processor
   */

    vd.ghost_get<type, rho, Pressure, Pressure_prev, velocity, velocity_prev, pressure_acc, viscous_acc>();

    auto part = vd.getDomainIterator();

    while (part.isNext())
    {
        auto a = part.get();

        // skip boundaryParticles
        if (vd.getProp<type>(a) == BOUNDARY)
        {
            ++part;
            continue;
        };

        // Get particle properties a
        Point<3, double> xa = vd.getPos(a);
        Point<3, double> va = vd.getProp<velocity>(a);

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);
        double Pa = vd.getProp<Pressure>(a);

        // Initialize quantities for summation in loop
        Point<3, double> gradientSum_1 = {0, 0, 0};
        Point<3, double> gradientSum_2 = {0, 0, 0};
        Point<3, double> drScaled = {0, 0, 0};

        double kernel = 0;
        double densitySumation = 0.0;
        double densityDiffusionTerm = 0.0;
        double pressureNeighbouring = 0.0;
        double sumSquaredGradient = 0.0;

        // get Neighborhood of particle
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            // ... q
            auto b = Np.get();

            // Get particle properties b

            Point<3, double> xb = vd.getPos(b);
            Point<3, double> vb = vd.getProp<velocity>(b);

            double Pb = vd.getProp<Pressure>(b);
            double rhob = vd.getProp<rho>(b);
            double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;

            // calculate quantities
            Point<3, double> dr = xa - xb;
            Point<3, double> v_rel = va - vb;
            Point<3, double> DW;

            double r2 = norm2(dr);
            double r = sqrt(r2);
            DWab(dr, DW, r, false);

            // If the particles interact ...
            if (r < smoothingRadius)
            {
                double wab = Wab(r);
                kernel += wab;
                // pressureLogger << "r  :         " << r << "\n";
                // pressureLogger << "H  :         " << H << "\n";
                // pressureLogger << "r/H  :         " << r/H << "\n";
                
                // pressureLogger << "wab  :         " << wab << "\n";
                // PressureCoefficient Summation Terms
                sumSquaredGradient += dotProduct(DW, DW);
                gradientSum_1 += DW;
                gradientSum_2 += DW;
            }
            ++Np;
        }
        // Pressure calculation
        maxShephard = std::max(maxShephard, kernel);
        double beta_candiate = dotProduct(gradientSum_1, gradientSum_2) + sumSquaredGradient;
        beta = (kernel > maxShephard) ? beta_candiate : beta;

        pressureLogger << "kernel  :         " << kernel << "\n";
        pressureLogger << "maxShephard  :         " << maxShephard << "\n";
        pressureLogger << "beta  :         " << beta << "\n";

        ++part;
    }
}

template <typename CellList>
inline void predictPressure(particles &vd, CellList &NN, double &max_visc, const double &dt, double &errorGlobal, double &beta)
{
    /*
     paralelization housekeeping
     currently irrelevant
     because testing with one processor
   */

    vd.ghost_get<type, rho, Pressure, Pressure_prev, velocity, velocity_prev, pressure_acc, viscous_acc>();

    auto part = vd.getDomainIterator();

    while (part.isNext())
    {
        auto a = part.get();

        // skip boundaryParticles
        if (vd.getProp<type>(a) == BOUNDARY)
        {
            ++part;
            continue;
        };

        // Get particle properties a
        Point<3, double> xa = vd.getPos(a);
        Point<3, double> va = vd.getProp<velocity>(a);

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);
        double Pa = vd.getProp<Pressure>(a);

        // Initialize quantities for summation in loop
        Point<3, double> gradientSum_1 = {0, 0, 0};
        Point<3, double> gradientSum_2 = {0, 0, 0};
        Point<3, double> drScaled = {0, 0, 0};

        double kernel = 0;
        double densitySumation = 0.0;
        double densityDiffusionTerm = 0.0;
        double pressureNeighbouring = 0.0;
        double sumSquaredGradient = 0.0;

        // get Neighborhood of particle
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            // ... q
            auto b = Np.get();

            // Get particle properties b

            Point<3, double> xb = vd.getPos(b);
            Point<3, double> vb = vd.getProp<velocity>(b);

            double Pb = vd.getProp<Pressure>(b);
            double rhob = vd.getProp<rho>(b);
            double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;

            // calculate quantities
            Point<3, double> dr = xa - xb;
            Point<3, double> v_rel = va - vb;
            Point<3, double> DW;

            double r2 = norm2(dr);
            double r = sqrt(r2);
            DWab(dr, DW, r, false);
            double wab = Wab(r);

            // If the particles interact ...
            if (r < smoothingRadius)
            {
                // PressureCoefficient Summation Terms
                sumSquaredGradient += massb / rhob * dotProduct(DW, DW);
                gradientSum_1 += massb / rhob * DW;
                gradientSum_2 += massb * DW;
            }
            ++Np;
        }
        // Pressure calculation

        // double beta = dotProduct(gradientSum_1, gradientSum_2) + massa * sumSquaredGradient;
        double delptaP = (rhoa * rhoa * vd.getProp<rho_err>(a)) / (massa * massa * dt * dt * beta);
        double predicted_Pressure = vd.getProp<Pressure>(a) + delptaP + Pbackground;

        if ((predicted_Pressure < 0.0))
        {
            pressureLogger << "predicted_Pressure :         " << predicted_Pressure << "\n";
            pressureLogger << "vd.getProp<Pressure>(a) :    " << vd.getProp<Pressure>(a) << "\n";
            pressureLogger << "beta :  " << beta << "\n";
            pressureLogger << "rhoError :  " << vd.getProp<rho_err>(a) << "\n";
            pressureLogger << "delptaP :  " << delptaP << "\n";
        }

        vd.getProp<Pressure>(a) = predicted_Pressure;
        // vd.getProp<Pressure>(a) = smoothedPressure;

        ++part;
    }
}
template <typename CellList>
inline void smoothPressure(particles &vd, CellList &NN, double &max_visc, const double &dt, double &errorGlobal)
{
    /*
     paralelization housekeeping
     currently irrelevant
     because testing with one processor
   */

    vd.ghost_get<type, rho, Pressure, Pressure_prev, velocity, velocity_prev, pressure_acc, viscous_acc>();

    auto part = vd.getDomainIterator();

    while (part.isNext())
    {
        auto a = part.get();

        // skip boundaryParticles
        if (vd.getProp<type>(a) == BOUNDARY)
        {
            ++part;
            continue;
        };

        // Get particle properties a
        Point<3, double> xa = vd.getPos(a);

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);
        double Pa = vd.getProp<Pressure>(a);

        double kernel = 0;
        double pressureNeighbouring = 0.0;

        // get Neighborhood of particle
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            // ... q
            auto b = Np.get();

            // Get particle properties b

            Point<3, double> xb = vd.getPos(b);

            double Pb = vd.getProp<Pressure>(b);
            double rhob = vd.getProp<rho>(b);
            double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;

            // calculate quantities
            Point<3, double> dr = xa - xb;

            double r2 = norm2(dr);
            double r = sqrt(r2);
            double wab = Wab(r);

            // If the particles interact ...
            if (r < smoothingRadius)
            {
                // Pressure Smoothing
                pressureNeighbouring += massb / rhob * vd.getProp<Pressure>(a) * wab;
                kernel += wab;
            }
            ++Np;
        }

        // Pressure calculation
        pressureNeighbouring *= (1.0 / kernel);

        // smoothed Pressure
        double smoothedPressure = vd.getProp<Pressure>(a) * intPConst + (1 - intPConst) * pressureNeighbouring;

        // TODO Pressure not clipped
        vd.getProp<Pressure>(a) = (kernel > EPS) ? smoothedPressure : Pbackground;
        // vd.getProp<Pressure>(a) = smoothedPressure;

        ++part;
    }
}

template <typename CellList>
inline void calc_Pressure_forces(particles &vd, CellList &NN, double &max_visc)
{
    vd.ghost_get<type, rho, Pressure, Pressure_prev, velocity, velocity_prev, pressure_acc, viscous_acc>();

    auto part = vd.getDomainIterator();

    while (part.isNext())
    {

        auto a = part.get();

        // skip boundary particles
        if (vd.getProp<type>(a) == BOUNDARY)
        {
            ++part;
            continue;
        };

        // Get particle properties a
        Point<3, double> xa = vd.getPos(a);

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);
        double Pa = vd.getProp<Pressure>(a);

        // get Neighborhood of particle
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        // For each neighborhood particle
        while (Np.isNext() == true)
        {

            auto b = Np.get();

            // Get particle propertiesb
            Point<3, double> xb = vd.getPos(b);

            double Pb = vd.getProp<Pressure>(b);
            double rhob = vd.getProp<rho>(b);
            double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;

            // calculate quantities
            Point<3, double> dr = xa - xb;
            double r2 = norm2(dr);
            double r = sqrt(r2);

            // If the particles interact ...
            if (r < smoothingRadius)
            {
                // DW = partial_W/partial_w_ab *e_ab --> Vector
                Point<3, double> DW;
                DWab(dr, DW, r, false);

                double factor = -massb * ((Pa + Pb) / (rhoa * rhob));

                vd.getProp<pressure_acc>(a)[0] += factor * DW.get(0);
                vd.getProp<pressure_acc>(a)[1] += factor * DW.get(1);
                vd.getProp<pressure_acc>(a)[2] += factor * DW.get(2);
            }
            ++Np;
        }
        ++part;
    }
}

template <typename CellList>
inline void extrapolateBoundaries(particles &vd, CellList &NN, double &max_visc)
{
    auto part = vd.getDomainIterator();

    vd.map();
    vd.ghost_get<type, rho, Pressure, Pressure_prev, velocity, velocity_prev, pressure_acc, viscous_acc>();
    vd.updateCellList(NN);

    while (part.isNext())
    {
        auto a = part.get();

        // skip Particle if it is Fluid
        if (vd.getProp<type>(a) == FLUID)
        {
            ++part;
            continue;
        }
        // Get particle properties a
        Point<3, double> xa = vd.getPos(a);
        Point<3, double> va = vd.getProp<velocity>(a);

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);
        double Pa = vd.getProp<Pressure>(a);

        // Initialize quantities for summation in loop

        Point<3, double> rho_dr_wf = {0, 0, 0};
        Point<3, double> v_wab_vectorial = {0, 0, 0};

        double kernel = 0;
        double rho_wab = 0;
        double p_wab = 0;
        double g_wab = 0;
        double v_wab = 0;

        // get Neighborhood of particle
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            auto b = Np.get();
            // Skip boundary particles only consider FLUID particles
            // for extrapolation of boundary condition
            if (vd.getProp<type>(b) == BOUNDARY)
            {
                ++Np;
                continue;
            };
            // Get particle properties b

            Point<3, double> xb = vd.getPos(b);
            Point<3, double> vb = vd.getProp<velocity>(b);

            double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
            double Pb = vd.getProp<Pressure>(b);
            double rhob = vd.getProp<rho>(b);

            // calculate quantities
            Point<3, double> dr = xa - xb;

            double r2 = norm2(dr);
            double r = sqrt(r2);
            double wab = Wab(r);

            // If the particles interact ...
            if (r < smoothingRadius)
            {
                kernel += wab;

                // pressure and bodyforce summation
                p_wab += Pb * wab;
                rho_dr_wf += wab * rhob * dr;

                // density sumation
                rho_wab += rhob * wab;

                // velocity sumation
                v_wab_vectorial += vb * wab;
            }
            ++Np;
        }

        // pressure extrapolation

        // scalar product g_ext and rho_dr_wf

        // hack that dotProduct()- function works:
        // qickAndDirty
        Point<3, double> BF{bodyforce[0], bodyforce[1], bodyforce[2]};

        g_wab = dotProduct(rho_dr_wf, BF);

        // calculate candidate pressure
        // double cand_pressure = (p_wab + g_wab) / (kernel + EPS);
        double cand_pressure = (p_wab) / (kernel + EPS);
        //  if kernel empty set pressure zero, otherwise cand pressure

        //  TODO Boundary Pressure & Density
        // vd.getProp<Pressure>(a) = (kernel > EPS) ? cand_pressure : 0.0;

        vd.getProp<Pressure>(a) = (kernel > EPS) ? cand_pressure : Pbackground;

        // calculate extrapolated densities
        double cand_density = rho_wab / (kernel + EPS);

        // if kernel is empty set density to reference density
        vd.getProp<rho>(a) = (kernel > EPS) ? cand_density : rho_zero;

        // set negative velocity
        vd.template getProp<velocity>(a)[0] = (kernel > EPS) ? -v_wab_vectorial[0] / (kernel + EPS) : 0;
        vd.template getProp<velocity>(a)[1] = (kernel > EPS) ? -v_wab_vectorial[1] / (kernel + EPS) : 0;
        vd.template getProp<velocity>(a)[2] = (kernel > EPS) ? -v_wab_vectorial[2] / (kernel + EPS) : 0;

        ++part;
    }
}

void euler_int_PCISPH(particles &vd, double dt)
{
    vd.ghost_get<type, rho, Pressure, Pressure_prev, velocity, velocity_prev, pressure_acc, viscous_acc>();

    // list of the particle to remove
    to_remove.clear();
    // particle iterator
    auto part = vd.getDomainIterator();
    double dt205 = dt * dt * 0.5;
    double dt2 = dt * 2.0;
    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();
        // if the particle is boundary
        if (vd.template getProp<type>(a) == BOUNDARY)
        {

            // Update rho
            double rhop = vd.template getProp<rho>(a);
            // Update only the density
            vd.template getProp<velocity>(a)[0] = 0.0;
            vd.template getProp<velocity>(a)[1] = 0.0;
            vd.template getProp<velocity>(a)[2] = 0.0;
            // vd.template getProp<Pressure>(a) = 0.0;

            vd.template getProp<rho_prev>(a) = rhop;
            ++part;
            continue;
        }
        //  TODO Nochmal nacschauen ob diese Integration auch Stroemer Verlet Schema entspricht

        // integrate velocity
        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + (vd.template getProp<pressure_acc>(a)[0] + vd.template getProp<viscous_acc>(a)[0] + bodyforce[0]) * dt;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + (vd.template getProp<pressure_acc>(a)[1] + vd.template getProp<viscous_acc>(a)[1] + bodyforce[1]) * dt;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + (vd.template getProp<pressure_acc>(a)[2] + vd.template getProp<viscous_acc>(a)[2] + bodyforce[2]) * dt;
        // integrate position
        vd.getPos(a)[0] = vd.template getProp<x_pre>(a)[0] + (vd.template getProp<velocity_prev>(a)[0] + vd.template getProp<velocity>(a)[0]) * dt * 0.5;
        vd.getPos(a)[1] = vd.template getProp<x_pre>(a)[1] + (vd.template getProp<velocity_prev>(a)[1] + vd.template getProp<velocity>(a)[1]) * dt * 0.5;
        vd.getPos(a)[2] = vd.template getProp<x_pre>(a)[2] + (vd.template getProp<velocity_prev>(a)[2] + vd.template getProp<velocity>(a)[2]) * dt * 0.5;

        // integrate density
        // TODO vd.template getProp<rho>(a) = (CandidateDensity > rho_zero) ? CandidateDensity : rho_zero;
        double CandidateDensity = vd.template getProp<rho_prev>(a) + dt * vd.template getProp<drho>(a);

        // With or without Clipping
        // vd.template getProp<rho>(a) = (CandidateDensity > rho_zero) ? CandidateDensity : rho_zero;
        // vd.getProp<Pressure>(a) = (vd.getProp<Pressure>(a) > 0.0) ? vd.getProp<Pressure>(a) : 0.0;

        // vd.template getProp<rho>(a) = CandidateDensity;

        // save velocity_prev and rho_prev
        vd.template getProp<velocity_prev>(a)[0] = vd.template getProp<velocity>(a)[0];
        vd.template getProp<velocity_prev>(a)[1] = vd.template getProp<velocity>(a)[1];
        vd.template getProp<velocity_prev>(a)[2] = vd.template getProp<velocity>(a)[2];
        vd.template getProp<rho_prev>(a) = vd.template getProp<rho>(a);
        vd.template getProp<x_pre>(a)[0] = vd.getPos(a)[0];
        vd.template getProp<x_pre>(a)[1] = vd.getPos(a)[1];
        vd.template getProp<x_pre>(a)[2] = vd.getPos(a)[2];
        vd.template getProp<Pressure_prev>(a) = vd.getProp<Pressure>(a);

        ++part;
    }
    // remove the particles
    vd.remove(to_remove, 0);
    // increment the iteration counter
    cnt++;
};