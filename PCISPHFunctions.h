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
        vd.template getProp<pressure_force>(a)[0] = 0.0;
        vd.template getProp<pressure_force>(a)[1] = 0.0;
        vd.template getProp<pressure_force>(a)[2] = 0.0;

        vd.template getProp<viscousFoces>(a)[0] = 0.0;
        vd.template getProp<viscousFoces>(a)[1] = 0.0;
        vd.template getProp<viscousFoces>(a)[2] = 0.0;

        vd.template getProp<Pressure>(a) = 0.0;

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
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                };

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

            vd.getProp<viscousFoces>(a)[0] = (kernel > 0.0) ? viscousFocesTermx : 0.0;
            vd.getProp<viscousFoces>(a)[1] = (kernel > 0.0) ? viscousFocesTermy : 0.0;
            vd.getProp<viscousFoces>(a)[2] = (kernel > 0.0) ? viscousFocesTermz : 0.0;
        }
        ++part;
    }
}


void position_and_velocity_prediction(particles &vd, double dt)
{
    to_remove.clear();
    auto part = vd.getDomainIterator();

    while (part.isNext())
    {
        // ... a
        auto a = part.get();

        if (vd.template getProp<type>(a) == BOUNDARY)
        {
            ++part;
            continue;
        }

        // integrate velocity
        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + (vd.template getProp<pressure_force>(a)[0] + vd.template getProp<viscousFoces>(a)[0] + bodyforce[0]) * dt;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + (vd.template getProp<pressure_force>(a)[1] + vd.template getProp<viscousFoces>(a)[1] + bodyforce[1]) * dt;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + (vd.template getProp<pressure_force>(a)[2] + vd.template getProp<viscousFoces>(a)[2] + bodyforce[2]) * dt;
        // integrate position
        vd.getPos(a)[0] = vd.template getProp<x_pre>(a)[0] + (vd.template getProp<velocity>(a)[0]) * dt;
        vd.getPos(a)[1] = vd.template getProp<x_pre>(a)[1] + (vd.template getProp<velocity>(a)[1]) * dt;
        vd.getPos(a)[2] = vd.template getProp<x_pre>(a)[2] + (vd.template getProp<velocity>(a)[2]) * dt;

        ++part;

        if (vd.getPos(a)[0] < Xmin || vd.getPos(a)[1] < Ymin || vd.getPos(a)[2] < Zmin ||
            vd.getPos(a)[0] > Xmax || vd.getPos(a)[1] > Ymax || vd.getPos(a)[2] > Zmax)
        {
            to_remove.add(a.getKey());
            std::cout << "added to remove" << std::endl;
        }
    }

    vd.remove(to_remove, 0);
    cnt++;
}


template <typename CellList>
inline void predictDensityAndUpdate(particles &vd, CellList &NN, double &max_visc, const double &dt, double &errorGlobal)
{
    vd.map();
    vd.ghost_get<>();
    vd.updateCellList(NN);

    auto part = vd.getDomainIterator();

    while (part.isNext())
    {
        // ... a
        auto a = part.get();

        Point<3, double> xa = vd.getPos(a);
        Point<3, double> va = vd.getProp<velocity>(a);

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);
        double Pa = vd.getProp<Pressure>(a);

        if (vd.getProp<type>(a) != FLUID)
        {
            ++part;
            continue;
        };

        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        double kernel = 0;
        vd.getProp<drho>(a) = 0.0;
        double densityTerm1 = 0.0;
        double densityTerm2 = 0.0;
        double pressureNeighbouring = 0.0;
        double sumSquaredGradient = 0.0;
        Point<3, double> gradientSum_1 = {0, 0, 0};
        Point<3, double> gradientSum_2 = {0, 0, 0};
        Point<3, double> drScaled = {0, 0, 0};

        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            // ... q
            auto b = Np.get();
            Point<3, double> xb = vd.getPos(b);
            if (a.getKey() == b)
            {
                ++Np;
                continue;
            };

            Point<3, double> vb = vd.getProp<velocity>(b);
            Point<3, double> dr = xa - xb;
            double Pb = vd.getProp<Pressure>(b);
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

                // Desnsity Predicition Terms
                double dot_vrel_DW = v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2);
                drScaled = 2 * massb * (vd.getProp<rho_prev>(a) + vd.getProp<rho_prev>(b)) / (vd.getProp<rho_prev>(b) * (r2 + Eta2)) * dr;

                densityTerm2 += drScaled.get(0) * DW.get(0) + drScaled.get(1) * DW.get(1) + drScaled.get(2) * DW.get(2);
                densityTerm1 += massb * dot_vrel_DW;

                // Desnsity Predicition Terms

                sumSquaredGradient += massb/rhob*(DW.get(0) * DW.get(0) + DW.get(1) * DW.get(1) + DW.get(2) * DW.get(2));
                gradientSum_1 += massb/rhob*DW;
                gradientSum_2 += massb*DW;

                pressureNeighbouring += massb / rhob * Pb * wab;
            }
            ++Np;
        }
        outFile << "densityTerm1 " << densityTerm1 << std::endl;
        outFile << "visco * H * cbar  " << visco * H * cbar * Eta2 << std::endl;
        outFile << "  densityTerm2 " << densityTerm2 << std::endl;

        vd.getProp<drho>(a) = densityTerm1;// + visco * H * cbar * densityTerm2;
        vd.getProp<rho>(a) = vd.getProp<rho>(a) + vd.getProp<drho>(a) * dt;

        double rho_error = vd.getProp<rho>(a) - rho_zero;
        // TODO  double rho_error = std::max(0,vd.getProp<rho>(a) - rho_zero); // try clipping error
        double rho_error_normalized = rho_error / rho_zero;

        errorGlobal = std::max(errorGlobal, rho_error_normalized);

        double gradientSum_squared = gradientSum_1.get(0) * gradientSum_2.get(0) + gradientSum_1.get(1) * gradientSum_2.get(1) + gradientSum_1.get(2) * gradientSum_2.get(2);
        double beta = gradientSum_squared + massa*sumSquaredGradient  ;

        double dt2 = dt * dt;
        double massa2 = massa * massa;
        double rho_zero2 = rho_zero * rho_zero;
        double K_i = rhoa / (dt2 * beta);

        double delptaP = K_i * rho_error;

        double predicted_Pressure = vd.getProp<Pressure>(a) + delptaP;

        double interpolP = predicted_Pressure * intPConst + (1 - intPConst) * pressureNeighbouring;
        // TODO fix Pressure Interpolation
        vd.getProp<Pressure>(a) = predicted_Pressure;

        ++part;
    }
}

template <typename CellList>
inline void calc_Pressure_forces(particles &vd, CellList &NN, double &max_visc)
{
    auto part = vd.getDomainIterator();

    vd.updateCellList(NN);

    while (part.isNext())
    {
        // ... a
        auto a = part.get();

        // Point<3, double> xa = vd.getPos(a);
        Point<3, double> xa = vd.getPos(a);

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);
        double Pa = vd.getProp<Pressure>(a);

        if (vd.getProp<type>(a) != FLUID)
        {
            ++part;
            continue;
        };

        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

        while (Np.isNext() == true)
        {
            // ... q
            auto b = Np.get();
            Point<3, double> xb = vd.getPos(b);

            if (a.getKey() == b)
            {
                ++Np;
                continue;
            };
            Point<3, double> dr = xa - xb;

            double Pb = vd.getProp<Pressure>(b);
            double rhob = vd.getProp<rho>(b);
            double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
            double r2 = norm2(dr);
            double r = sqrt(r2);

            // If the particles interact ...
            if (r < smoothingRadius)
            {

                Point<3, double> DW;

                DWab(dr, DW, r, false);

                double factor = -massb * ((Pa + Pb) / (rhoa * rhob));

                vd.getProp<pressure_force>(a)[0] += factor * DW.get(0);
                vd.getProp<pressure_force>(a)[1] += factor * DW.get(1);
                vd.getProp<pressure_force>(a)[2] += factor * DW.get(2);
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
    vd.updateCellList(NN);

    while (part.isNext())
    {
        auto a = part.get();

        Point<3, double> xa = vd.getPos(a);
        Point<3, double> va = vd.getProp<velocity>(a);

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);
        double Pa = vd.getProp<Pressure>(a);

        if (vd.getProp<type>(a) != BOUNDARY)
        {
            ++part;
            continue;
        }
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        double kernel = 0;
        double rho_wab = 0;
        double p_wab = 0;
        double g_wab = 0;
        double v_wab = 0;

        Point<3, double> g_wab_vectorial = {0, 0, 0};
        Point<3, double> v_wab_vectorial = {0, 0, 0};
        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            // ... q
            auto b = Np.get();

            Point<3, double> xb = vd.getPos(b);
            Point<3, double> vb = vd.getProp<velocity>(b);

            double Pb = vd.getProp<Pressure>(b);

            if (a.getKey() == b || vd.getProp<type>(b) != FLUID)
            {
                ++Np;
                continue;
            };

            // get the mass of the particle
            double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
            double rhob = vd.getProp<rho>(b);
            Point<3, double> dr = xa - xb;

            double r2 = norm2(dr);
            double r = sqrt(r2);

            // If the particles interact ...
            if (r < smoothingRadius)
            {

                double wab = Wab(r);
                kernel += wab;
                rho_wab += rhob * wab;
                p_wab += Pb * wab;
                v_wab_vectorial += vb * wab;
                g_wab_vectorial += wab * rhob * dr;
                // Pressure and density calclation Boundary Particle
            }

            ++Np;
        }
        g_wab = g_wab_vectorial.get(0) * bodyforce[0] + g_wab_vectorial.get(1) * bodyforce[1] + g_wab_vectorial.get(2) * bodyforce[2];

        double cand_pressure = (p_wab + g_wab) / kernel;
        // double cand_pressure = (p_wab) / kernel;
        double cand_density = rho_wab / kernel;

        vd.getProp<rho>(a) = (kernel > 0.0) ? cand_density : rho_zero;
        vd.getProp<Pressure>(a) = (kernel > 0.0) ? cand_pressure : 0.0;
        // vd.getProp<Pressure>(a) = cand_pressure ;

        vd.template getProp<velocity>(a)[0] = (kernel > 0.0) ? -v_wab_vectorial[0] / kernel : 0;
        vd.template getProp<velocity>(a)[1] = (kernel > 0.0) ? -v_wab_vectorial[1] / kernel : 0;
        vd.template getProp<velocity>(a)[2] = (kernel > 0.0) ? -v_wab_vectorial[2] / kernel : 0;
        ++part;
    }
}

void stroemer_verlet_int(particles &vd, double dt)
{
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

        // integrate velocity
        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + (vd.template getProp<pressure_force>(a)[0] + vd.template getProp<viscousFoces>(a)[0] + bodyforce[0]) * dt;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + (vd.template getProp<pressure_force>(a)[1] + vd.template getProp<viscousFoces>(a)[1] + bodyforce[1]) * dt;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + (vd.template getProp<pressure_force>(a)[2] + vd.template getProp<viscousFoces>(a)[2] + bodyforce[2]) * dt;
        // integrate position
        vd.getPos(a)[0] = vd.template getProp<x_pre>(a)[0] + (vd.template getProp<velocity_prev>(a)[0] + vd.template getProp<velocity>(a)[0]) * dt * 0.5;
        vd.getPos(a)[1] = vd.template getProp<x_pre>(a)[1] + (vd.template getProp<velocity_prev>(a)[1] + vd.template getProp<velocity>(a)[1]) * dt * 0.5;
        vd.getPos(a)[2] = vd.template getProp<x_pre>(a)[2] + (vd.template getProp<velocity_prev>(a)[2] + vd.template getProp<velocity>(a)[2]) * dt * 0.5;

        // integrate density
        double CandidateDensity = vd.template getProp<rho_prev>(a) + dt * vd.template getProp<drho>(a);
        // vd.template getProp<rho>(a) = (CandidateDensity > rho_zero) ? CandidateDensity : rho_zero;
        vd.template getProp<rho>(a) = CandidateDensity;

        // save velocity_prev and rho_prev
        vd.template getProp<velocity_prev>(a)[0] = vd.template getProp<velocity>(a)[0];
        vd.template getProp<velocity_prev>(a)[1] = vd.template getProp<velocity>(a)[1];
        vd.template getProp<velocity_prev>(a)[2] = vd.template getProp<velocity>(a)[2];
        vd.template getProp<rho_prev>(a) = vd.template getProp<rho>(a);
        vd.template getProp<x_pre>(a)[0] = vd.getPos(a)[0];
        vd.template getProp<x_pre>(a)[1] = vd.getPos(a)[1];
        vd.template getProp<x_pre>(a)[2] = vd.getPos(a)[2];

        ++part;
    }
    // remove the particles
    vd.remove(to_remove, 0);
    // increment the iteration counter
    cnt++;
};