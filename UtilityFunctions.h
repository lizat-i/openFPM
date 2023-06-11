#pragma once
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include "inputs.h"
#include "vanillaFunctions.h"

particles InitializeDomain_vanillaSPH()
{
    // It contain for each time-step the value detected by the probes
    openfpm::vector<openfpm::vector<double>> press_t;

    // Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
    Box<3, double> domain({-0.05, -0.05, -0.05}, {1.7010, 0.7065, 0.5025});
    size_t sz[3] = {207, 90, 66};
    // Fill W_dap
    W_dap = 1.0 / Wab(H / 1.5);
    // Here we define the boundary conditions of our problem
    size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
    // extended boundary around the domain, and the processor domain
    Ghost<3, double> g(2 * H);

    particles vd(0, domain, bc, g, DEC_GRAN(512));
    // You can ignore all these dp/2.0 is a trick to reach the same initialization
    // of Dual-SPH that use a different criteria to draw particles

    Box<3, double> fluid_box({0, 0, 0}, {1, 1, 1});
    // return an iterator to the fluid particles to add to vd
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);
    // here we fill some of the constants needed by the simulation
    max_fluid_height = fluid_it.getBoxMargins().getHigh(2);
    h_swl = fluid_it.getBoxMargins().getHigh(2) - fluid_it.getBoxMargins().getLow(2);

    B = (coeff_sound) * (coeff_sound)*gravity * h_swl * rho_zero / gamma_;
    // cbar = coeff_sound * sqrt(gravity * h_swl);
    // for each particle inside the fluid box ...
    while (fluid_it.isNext())
    {
        // ... add a particle ...
        vd.add();
        // ... and set it position ...
        vd.getLastPos()[0] = fluid_it.get().get(0);
        vd.getLastPos()[1] = fluid_it.get().get(1);
        vd.getLastPos()[2] = fluid_it.get().get(2);
        // and its type.
        vd.template getLastProp<type>() = FLUID;
        // We also initialize the density of the particle and the hydro-static pressure given by
        //
        // rho_zero*g*h = P
        //
        // rho_p = (P/B + 1)^(1/Gamma) * rho_zero
        //
        vd.template getLastProp<Pressure>() = rho_zero * gravity * (max_fluid_height - fluid_it.get().get(2));
        vd.template getLastProp<rho>() = pow(vd.template getLastProp<Pressure>() / B + 1, 1.0 / gamma_) * rho_zero;
        vd.template getLastProp<rho_prev>() = vd.template getLastProp<rho>();
        vd.template getLastProp<velocity>()[0] = 0.0;
        vd.template getLastProp<velocity>()[1] = 0.0;
        vd.template getLastProp<velocity>()[2] = 0.0;
        vd.template getLastProp<velocity_prev>()[0] = 0.0;
        vd.template getLastProp<velocity_prev>()[1] = 0.0;
        vd.template getLastProp<velocity_prev>()[2] = 0.0;
        // next fluid particle
        ++fluid_it;
    }
    // Recipient
    Box<3, double> recipient1({0.0, 0.0, 0.0}, {1.6 + dp / 2.0, 0.67 + dp / 2.0, 0.4 + dp / 2.0});
    Box<3, double> recipient2({dp, dp, dp}, {1.6 - dp / 2.0, 0.67 - dp / 2.0, 0.4 + dp / 2.0});
    Box<3, double> obstacle1({0.9, 0.24 - dp / 2.0, 0.0}, {1.02 + dp / 2.0, 0.36, 0.45 + dp / 2.0});
    Box<3, double> obstacle2({0.9 + dp, 0.24 + dp / 2.0, 0.0}, {1.02 - dp / 2.0, 0.36 - dp, 0.45 - dp / 2.0});
    Box<3, double> obstacle3({0.9 + dp, 0.24, 0.0}, {1.02, 0.36, 0.45});
    openfpm::vector<Box<3, double>> holes;
    holes.add(recipient2);
    holes.add(obstacle1);
    auto bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient1);
    while (bound_box.isNext())
    {
        vd.add();
        vd.getLastPos()[0] = bound_box.get().get(0);
        vd.getLastPos()[1] = bound_box.get().get(1);
        vd.getLastPos()[2] = bound_box.get().get(2);
        vd.template getLastProp<type>() = BOUNDARY;
        vd.template getLastProp<rho>() = rho_zero;
        vd.template getLastProp<rho_prev>() = rho_zero;
        vd.template getLastProp<velocity>()[0] = 0.0;
        vd.template getLastProp<velocity>()[1] = 0.0;
        vd.template getLastProp<velocity>()[2] = 0.0;
        vd.template getLastProp<velocity_prev>()[0] = 0.0;
        vd.template getLastProp<velocity_prev>()[1] = 0.0;
        vd.template getLastProp<velocity_prev>()[2] = 0.0;
        ++bound_box;
    }
    auto obstacle_box = DrawParticles::DrawSkin(vd, sz, domain, obstacle2, obstacle1);
    while (obstacle_box.isNext())
    {
        vd.add();
        vd.getLastPos()[0] = obstacle_box.get().get(0);
        vd.getLastPos()[1] = obstacle_box.get().get(1);
        vd.getLastPos()[2] = obstacle_box.get().get(2);
        vd.template getLastProp<type>() = BOUNDARY;
        vd.template getLastProp<rho>() = rho_zero;
        vd.template getLastProp<rho_prev>() = rho_zero;
        vd.template getLastProp<velocity>()[0] = 0.0;
        vd.template getLastProp<velocity>()[1] = 0.0;
        vd.template getLastProp<velocity>()[2] = 0.0;
        vd.template getLastProp<velocity_prev>()[0] = 0.0;
        vd.template getLastProp<velocity_prev>()[1] = 0.0;
        vd.template getLastProp<velocity_prev>()[2] = 0.0;
        ++obstacle_box;
    }

    return vd;
}

particles InitializeDomain_PoiseuilleFlow3d(){
    // Initialize domain
    Box<3, double> domain({0, -4*dp, 0}, {1, 1+4*dp, 1});
    size_t sz[3] = {NumberOfParticles_x, NumberOfParticles_y, NumberOfParticles_z};
    size_t bc[3] = {PERIODIC, NON_PERIODIC, PERIODIC};
    Ghost<3, double> g(2 * H);
    particles vd(0, domain, bc, g, DEC_GRAN(512));



    Box<3, double> fluid_box({0, 0, 0}, {1, 1, 1});
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

    while (fluid_it.isNext())
    {

        vd.add();
        vd.getLastPos()[0] = fluid_it.get().get(0);
        vd.getLastPos()[1] = fluid_it.get().get(1);
        vd.getLastPos()[2] = fluid_it.get().get(2);

        vd.template getLastProp<type>() = FLUID;
        vd.template getLastProp<Pressure>() = 0;

        vd.template getLastProp<drho>() = 0;
        vd.template getLastProp<rho>() = rho_zero;
        vd.template getLastProp<rho_prev>() = rho_zero;

        vd.template getLastProp<velocity>()[0] = 0.0;
        vd.template getLastProp<velocity>()[1] = 0.0;
        vd.template getLastProp<velocity>()[2] = 0.0;

        vd.template getLastProp<force>()[0] = 0.0;
        vd.template getLastProp<force>()[1] = 0.0;
        vd.template getLastProp<force>()[2] = 0.0;

        vd.template getLastProp<velocity_prev>()[0] = 0.0;
        vd.template getLastProp<velocity_prev>()[1] = 0.0;
        vd.template getLastProp<velocity_prev>()[2] = 0.0;
        ++fluid_it;
    }

    // Recipient
    Box<3, double> yminplate({0.0, 0.0 - 4 * dp, 0.0}, {1.0, 0.0, 1.0});
    Box<3, double> ymaxplate({0.0, 1.0, 0.0}, {1.0, 1.0 + 4 * dp, 1.0});

    openfpm::vector<Box<3, double>> BoundariesBoxes;
    BoundariesBoxes.add(yminplate);
    BoundariesBoxes.add(ymaxplate);

    for (size_t i = 0; i < 2; ++i)
    {
        std::cout << "for box number " << i << '\n';
        Box<3, double> box = BoundariesBoxes.get(i);
        auto toFill = DrawParticles::DrawBox(vd, sz, domain, box);

        while (toFill.isNext())
        {

            vd.add();
            vd.getLastPos()[0] = toFill.get().get(0);
            vd.getLastPos()[1] = toFill.get().get(1);
            vd.getLastPos()[2] = toFill.get().get(2);

            vd.template getLastProp<type>() = BOUNDARY;
            vd.template getLastProp<Pressure>() = 0;

            vd.template getLastProp<drho>() = 0;
            vd.template getLastProp<rho>() = rho_zero;
            vd.template getLastProp<rho_prev>() = rho_zero;

            vd.template getLastProp<velocity>()[0] = 0.0;
            vd.template getLastProp<velocity>()[1] = 0.0;
            vd.template getLastProp<velocity>()[2] = 0.0;

            vd.template getLastProp<force>()[0] = 0.0;
            vd.template getLastProp<force>()[1] = 0.0;
            vd.template getLastProp<force>()[2] = 0.0;

            vd.template getLastProp<velocity_prev>()[0] = 0.0;
            vd.template getLastProp<velocity_prev>()[1] = 0.0;
            vd.template getLastProp<velocity_prev>()[2] = 0.0;
            ++toFill;
        }
    }
 
    return vd;
}