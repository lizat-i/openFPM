#pragma once
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include "inputs.h"
#include "vanillaFunctions.h"

particles InitializeDomain_vanillaSPH()
{

    Box<3, double> domain({-0.05, -0.05, -0.05}, {1.7010, 0.7065, 0.5025});
    size_t sz[3] = {207, 90, 66};
    W_dap = 1.0 / Wab(H / 1.5);
    size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
    Ghost<3, double> g(2 * H);
    particles vd(0, domain, bc, g, DEC_GRAN(512));
    Box<3, double> fluid_box({dp / 2.0, dp / 2.0, dp / 2.0}, {0.4 + dp / 2.0, 0.67 - dp / 2.0, 0.3 + dp / 2.0});
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

    max_fluid_height = fluid_it.getBoxMargins().getHigh(2);
    h_swl = fluid_it.getBoxMargins().getHigh(2) - fluid_it.getBoxMargins().getLow(2);
    B = (coeff_sound) * (coeff_sound)*gravity * h_swl * rho_zero / gamma_;
 
    int NrOfFluidParticles = 0;
    while (fluid_it.isNext())
    {
        vd.add();
        vd.getLastPos()[0] = fluid_it.get().get(0);
        vd.getLastPos()[1] = fluid_it.get().get(1);
        vd.getLastPos()[2] = fluid_it.get().get(2);

        vd.template getLastProp<x_pre>()[0] = vd.getLastPos()[0];
        vd.template getLastProp<x_pre>()[1] = vd.getLastPos()[1];
        vd.template getLastProp<x_pre>()[2] = vd.getLastPos()[2];

        vd.template getLastProp<type>() = FLUID;
        vd.template getLastProp<Pressure>() = rho_zero * gravity * (max_fluid_height - fluid_it.get().get(2));
        vd.template getLastProp<rho>() = pow(vd.template getLastProp<Pressure>() / B + 1, 1.0 / gamma_) * rho_zero;
        vd.template getLastProp<rho_prev>() = vd.template getLastProp<rho>();
        vd.template getLastProp<velocity>()[0] = 0.0;
        vd.template getLastProp<velocity>()[1] = 0.0;
        vd.template getLastProp<velocity>()[2] = 0.0;
        vd.template getLastProp<velocity_prev>()[0] = 0.0;
        vd.template getLastProp<velocity_prev>()[1] = 0.0;
        vd.template getLastProp<velocity_prev>()[2] = 0.0;
        ++fluid_it;
        ++NrOfFluidParticles;
    }

    // Build domain with 6 boxes.
    // first groundplate
    Box<3, double> groundplate({0.0 - 3 * dp, 0.0 - 3 * dp, 0.0 - 3.0 * dp}, {1.6 + dp * 3.0, 0.67 + dp * 3.0, 0.0 + dp / 2});
    // xmin xmax plates
    Box<3, double> xmin({0.0 - 3 * dp, 0.0 - 3 * dp, 0.0}, {0.0 + dp / 2, 0.67 + dp * 3.0, 0.45 + dp / 2});
    Box<3, double> xmax({1.6 - dp / 2, 0.0 - 3 * dp, 0.0}, {1.6 + dp * 3.0, 0.67 + dp * 3.0, 0.45 + dp / 2});
    // ymin ymax plates
    Box<3, double> ymin({0.0 - 3 * dp, 0.0 - 3 * dp, 0.0}, {1.6 + dp * 3.0 + dp / 2, 0.00, 0.4 + dp / 2});
    Box<3, double> ymax({0.0 - 3 * dp, 0.67 + dp / 2, 0.0}, {1.6 + dp * 3.0, 0.67 + dp * 3.0, 0.4 + dp / 2});
    // obstacle box
    Box<3, double> obstacle1({0.9, 0.24 - dp / 2.0, 0.0}, {1.02 + dp / 2.0, 0.36, 0.45 + dp / 2});

    openfpm::vector<Box<3, double>> obstacle_and_bound_box;

    obstacle_and_bound_box.add(groundplate);
    obstacle_and_bound_box.add(xmin);
    obstacle_and_bound_box.add(xmax);
    obstacle_and_bound_box.add(ymin);
    obstacle_and_bound_box.add(ymax);
    obstacle_and_bound_box.add(obstacle1);
    std::cout << "draw box" << '\n';
    for (size_t i = 0; i < 6; ++i)
    {
        std::cout << "for box number " << i << '\n';
        Box<3, double> box = obstacle_and_bound_box.get(i);
        auto toFill = DrawParticles::DrawBox(vd, sz, domain, box);

        while (toFill.isNext())
        {
            vd.add();
            vd.getLastPos()[0] = toFill.get().get(0);
            vd.getLastPos()[1] = toFill.get().get(1);
            vd.getLastPos()[2] = toFill.get().get(2);

            vd.template getLastProp<x_pre>()[0] = vd.getLastPos()[0];
            vd.template getLastProp<x_pre>()[1] = vd.getLastPos()[1];
            vd.template getLastProp<x_pre>()[2] = vd.getLastPos()[2];

            vd.template getLastProp<type>() = BOUNDARY;
            vd.template getLastProp<rho>() = rho_zero;
            vd.template getLastProp<rho_prev>() = rho_zero;
            vd.template getLastProp<velocity>()[0] = 0.0;
            vd.template getLastProp<velocity>()[1] = 0.0;
            vd.template getLastProp<velocity>()[2] = 0.0;
            vd.template getLastProp<velocity_prev>()[0] = 0.0;
            vd.template getLastProp<velocity_prev>()[1] = 0.0;
            vd.template getLastProp<velocity_prev>()[2] = 0.0;
            ++toFill;
        }
    }
}

particles InitializeDomain_PoiseuilleFlow3d()
{
    // Initialize domain
    Box<3, double> domain({0, -4 * dp, -4 * dp}, {1, 1 + 4 * dp, 4 * dp});
    size_t sz[3] = {NumberOfParticles_x, NumberOfParticles_y, NumberOfParticles_z};
    size_t bc[3] = {PERIODIC, NON_PERIODIC, PERIODIC};
    Ghost<3, double> g(2 * H);
    particles vd(0, domain, bc, g, DEC_GRAN(512));
    W_dap = 1.0 / Wab(H / 1.5);
    Box<3, double> fluid_box({0, 0, -4 * dp}, {1, 1, 4 * dp});
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
    Box<3, double> yminplate({0.0, 0.0 - 4 * dp, -4 * dp}, {1.0, 0.0, 4 * dp});
    Box<3, double> ymaxplate({0.0, 1.0, -4 * dp}, {1.0, 1.0 + 4 * dp, 4 * dp});

    openfpm::vector<Box<3, double>> BoundariesBoxes;
    BoundariesBoxes.add(yminplate);
    BoundariesBoxes.add(ymaxplate);

    for (size_t i = 0; i < 2; ++i)
    {
        std::cout << "NumberOfParticles_z " << NumberOfParticles_z << '\n';
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