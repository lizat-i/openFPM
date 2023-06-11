#pragma once
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include "inputs.h"
#include "vanillaFunctions.h"
#include "UtilityFunctions.h"

void peng_int(particles &vd, double dt)
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
            // vd.template getProp<velocity>(a)[0] = 0.0;
            // vd.template getProp<velocity>(a)[1] = 0.0;
            // vd.template getProp<velocity>(a)[2] = 0.0;
            // double rhonew = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);
            // vd.template getProp<rho>(a) = (rhonew < rho_zero) ? rho_zero : rhonew;
            vd.template getProp<rho_prev>(a) = rhop;
            ++part;
            continue;
        }
        //-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
        double dx = vd.template getProp<velocity>(a)[0] * dt + (vd.template getProp<force>(a)[0] + vd.template getProp<viscousFoces>(a)[0]) * dt205;
        double dy = vd.template getProp<velocity>(a)[1] * dt + (vd.template getProp<force>(a)[1] + vd.template getProp<viscousFoces>(a)[1]) * dt205;
        double dz = vd.template getProp<velocity>(a)[2] * dt + (vd.template getProp<force>(a)[2] + vd.template getProp<viscousFoces>(a)[2]) * dt205;
        vd.getPos(a)[0] += dx;
        vd.getPos(a)[1] += dy;
        vd.getPos(a)[2] += dz;
        double velX = vd.template getProp<velocity>(a)[0];
        double velY = vd.template getProp<velocity>(a)[1];
        double velZ = vd.template getProp<velocity>(a)[2];
        double rhop = vd.template getProp<rho>(a);
        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity>(a)[0] + (vd.template getProp<force>(a)[0] + vd.template getProp<viscousFoces>(a)[0]) * dt;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity>(a)[1] + (vd.template getProp<force>(a)[1] + vd.template getProp<viscousFoces>(a)[1]) * dt;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity>(a)[2] + (vd.template getProp<force>(a)[2] + vd.template getProp<viscousFoces>(a)[2]) * dt;
        vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);
        // Check if the particle go out of range in space and in density
        if (vd.template getProp<rho>(a) < RhoMin || vd.template getProp<rho>(a) > RhoMax)
        {
            // to_remove.add(a.getKey());
        }
        vd.template getProp<velocity_prev>(a)[0] = velX;
        vd.template getProp<velocity_prev>(a)[1] = velY;
        vd.template getProp<velocity_prev>(a)[2] = velZ;
        vd.template getProp<rho_prev>(a) = rhop;
        ++part;
    }
    // remove the particles
    vd.remove(to_remove, 0);
    // increment the iteration counter
    cnt++;
}


template <typename CellList>
inline void calc_forces_boundary(particles &vd, CellList &NN, double &max_visc, auto &a, Point<3, double> &xa, Point<3, double> &va, double &massa, double &rhoa, double &Pa)
{
    vd.template getProp<force>(a)[0] = 0;
    vd.template getProp<force>(a)[1] = 0;
    vd.template getProp<force>(a)[2] = 0;

    vd.template getProp<viscousFoces>(a)[0] = 0;
    vd.template getProp<viscousFoces>(a)[1] = 0;
    vd.template getProp<viscousFoces>(a)[2] = 0;

    auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

    double kernel = 0;
    double rho_wab = 0;
    double p_wab = 0;
    double g_wab = 0;
    double v_wab = 0;

    Point<3, double> g_wab_vec = {0, 0, 0};
    Point<3, double> v_wab_vec = {0, 0, 0};

    while (Np.isNext() == true)
    {
        auto b = Np.get();

        if (a.getKey() == b || vd.getProp<type>(b) == BOUNDARY)
        {
            ++Np;
            continue;
        };

        Point<3, double> xb = vd.getPos(b);
        Point<3, double> dr = xa - xb;
        double r2 = norm2(dr);

        if (r2 < 4.0 * H * H)
        {

            double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
            double Pb = vd.getProp<Pressure>(b);
            double rhob = vd.getProp<rho>(b);
            double r = sqrt(r2);
            double wab = Wab(r);

            Point<3, double> vb = vd.getProp<velocity>(b);
            Point<3, double> dv = va - vb;
            Point<3, double> DW;
            double partialDWAB = 0;
            DWab(dr, DW, r, false, partialDWAB);

            const double dot = dr.get(0) * dv.get(0) + dr.get(1) * dv.get(1) + dr.get(2) * dv.get(2);

            const double dot_rr2 = dot / (r2 + Eta2);
            max_visc = std::max(dot_rr2, max_visc);

            kernel += wab;
            rho_wab += rhob * wab;
            p_wab += Pb * wab;
            g_wab_vec += wab * rhob  * dr ;
            v_wab_vec += vb * wab;

            vd.getProp<drho>(a) += massb * (dv.get(0) * DW.get(0) + dv.get(1) * DW.get(1) + dv.get(2) * DW.get(2));
        }
        ++Np;
    }
    g_wab = ( bodyforce[0] * g_wab_vec.get(0) +  bodyforce[1]  * g_wab_vec.get(1) + bodyforce[2] * g_wab_vec.get(2));

    double candidate_density = rho_wab / kernel;
    double candidate_Pressure = (p_wab + g_wab) / kernel;

    double candidate_velocity_x = -(v_wab_vec.get(0)) * (1 / kernel);
    double candidate_velocity_y = -(v_wab_vec.get(1)) * (1 / kernel);
    double candidate_velocity_z = -(v_wab_vec.get(2)) * (1 / kernel);

    vd.getProp<velocity>(a)[0] = (kernel > eps) ? candidate_velocity_x : 0.0;
    vd.getProp<velocity>(a)[1] = (kernel > eps) ? candidate_velocity_y : 0.0;
    vd.getProp<velocity>(a)[2] = (kernel > eps) ? candidate_velocity_z : 0.0;

    // TODO Boundary Clipping of Pressure and Density

    // with boundary pressure and density clipping

    vd.getProp<rho>(a) = (candidate_density > rho_zero) ? candidate_density : rho_zero;
    vd.getProp<Pressure>(a) = (candidate_Pressure > eps) ? candidate_Pressure : 0.0;

    // without boundary pressure and density clipping

    // vd.getProp<rho>(b) = candidate_density          ;
    // vd.getProp<Pressure>(b) = candidate_Pressure    ;
}

template <typename CellList>
inline void calc_forces_fluid(particles &vd, CellList &NN, double &max_visc, auto &a, Point<3, double> &xa, Point<3, double> &va, double &massa, double &rhoa, double &Pa)
{

    vd.template getProp<force>(a)[0] = bodyforce[0];
    vd.template getProp<force>(a)[1] = bodyforce[1];
    vd.template getProp<force>(a)[2] = bodyforce[2];

    vd.template getProp<viscousFoces>(a)[0] = 0;
    vd.template getProp<viscousFoces>(a)[1] = 0;
    vd.template getProp<viscousFoces>(a)[2] = 0;

    double viscousFocesTermx = 0;
    double viscousFocesTermy = 0;
    double viscousFocesTermz = 0;

    auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

    while (Np.isNext() == true)
    {

        auto b = Np.get();

        if (a.getKey() == b)
        {
            ++Np;
            continue;
        };

        Point<3, double> xb = vd.getPos(b);
        Point<3, double> vb = vd.getProp<velocity>(b);
        Point<3, double> dr = xa - xb;

        double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
        double Pb = vd.getProp<Pressure>(b);
        double rhob = vd.getProp<rho>(b);
        double r2 = norm2(dr);
        double Va = massa / rhoa;
        double Vb = massb / rhob;

        if (r2 < 4.0 * H * H)
        {
            double r = sqrt(r2);
            Point<3, double> v_rel = va - vb;
            Point<3, double> DW;
            double partialDWAB = 0;
            DWab(dr, DW, r, false, partialDWAB);

            double factor = -massb * ((vd.getProp<Pressure>(a) + vd.getProp<Pressure>(b)) / (rhoa * rhob) + Tensile(r, rhoa, rhob, Pa, Pb) ) ; //+ Pi(dr, r2, v_rel, rhoa, rhob, massb, max_visc));

            vd.getProp<force>(a)[0] += factor * DW.get(0);
            vd.getProp<force>(a)[1] += factor * DW.get(1);
            vd.getProp<force>(a)[2] += factor * DW.get(2);

            double viscousForcesFactor = -((2*kinematic_viscosity*kinematic_viscosity)/(kinematic_viscosity + kinematic_viscosity)) * (Va * Va + Vb * Vb)*(partialDWAB / r);

            viscousFocesTermx += viscousForcesFactor * v_rel.get(0);
            viscousFocesTermy += viscousForcesFactor * v_rel.get(1);
            viscousFocesTermz += viscousForcesFactor * v_rel.get(2);

            vd.getProp<drho>(a) += massb * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
        }
        ++Np;
    }
    vd.getProp<viscousFoces>(a)[0] = viscousFocesTermx/massa;
    vd.getProp<viscousFoces>(a)[1] = viscousFocesTermy/massa;
    vd.getProp<viscousFoces>(a)[2] = viscousFocesTermz/massa;
}

template <typename CellList>
inline void calc_forces_peng(particles &vd, CellList &NN, double &max_visc)
{

    auto part = vd.getDomainIterator();
    // Update the cell-list
    vd.updateCellList(NN);
    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();

        Point<3, double> xa = vd.getPos(a);
        Point<3, double> va = vd.getProp<velocity>(a);

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);
        double Pa = vd.getProp<Pressure>(a);

        vd.template getProp<drho>(a) = 0.0;

        if (vd.getProp<type>(a) != FLUID)
        {
            calc_forces_boundary(vd, NN, max_visc, a, xa, va, massa, rhoa, Pa);
        }
        else
        {
            calc_forces_fluid(vd, NN, max_visc, a, xa, va, massa, rhoa, Pa);
        }
        ++part;
    }
}
