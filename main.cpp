// #define SE_CLASS1
// #define STOP_ON_ERROR
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include <limits>

#include "main.h"

int main(int argc, char *argv[])
{
    std::cout << "started programm" << '\n';

    openfpm_init(&argc, &argv);
    openfpm::vector<openfpm::vector<double>> press_t;

    double Xmin = 0;
    double Ymin = 0 - dp * 3;
    double Zmin = 0;
    double Xmax = 1;
    double Ymax = 1 + dp * 3;
    double Zmax = 1;

    double L_x = Xmax - Xmin;
    double L_y = Ymax - Ymin;
    double L_z = Zmax - Zmin;

    double characteristic_Length = 1;

    size_t NumberOfParticles_x = ceil(L_x / dp);
    size_t NumberOfParticles_y = ceil(L_y / dp);
    size_t NumberOfParticles_z = ceil(L_z / dp);

    Box<3, double> domain({Xmin, Ymin, Zmin}, {Xmax, Ymax, Zmax});
    size_t sz[3] = {NumberOfParticles_x, NumberOfParticles_y, NumberOfParticles_z};
    W_dap = 1.0 / Wab(H / 1.5);

    size_t bc[3] = {PERIODIC, NON_PERIODIC, PERIODIC};
    Ghost<3, double> g(2 * H);
    particles vd(0, domain, bc, g, DEC_GRAN(512));
    Box<3, double> fluid_box({0, 0, 0}, {1.0, 1.0, 1.0});
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

    double cbar_veloc = (initialVelocity * initialVelocity) / densityVariation;
    double cbar_visco = dynamic_viscosity * initialVelocity / (densityVariation * characteristic_Length);
    double cbar_accel = gravity * characteristic_Length / densityVariation;
    cbar = std::max(cbar_veloc, cbar_visco);
    LOGDouble("cbar_veloc", cbar_veloc);
    LOGDouble("cbar_visco", cbar_visco);
    LOGDouble("cbar", cbar);

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
        vd.template getLastProp<Pressure>() = 0;
        vd.template getLastProp<rho>() = rho_zero;
        vd.template getLastProp<rho_prev>() = vd.template getLastProp<rho>();
        vd.template getLastProp<velocity>()[0] = 0.0;
        vd.template getLastProp<velocity>()[1] = 0.0;
        vd.template getLastProp<velocity>()[2] = 0.0;
        vd.template getLastProp<velocity_prev>()[0] = 0.0;
        vd.template getLastProp<velocity_prev>()[1] = 0.0;
        vd.template getLastProp<velocity_prev>()[2] = 0.0;
        ++fluid_it;
    }

    // Build domain with 6 boxes.
    // first groundplate

    Box<3, double> yminplate({0.0 - 3 * dp, 0.0 - 3 * dp, 0.0 - 3 * dp}, {1.0 + 3 * dp, 0.0, 1.0 + 3 * dp});
    Box<3, double> ymaxplate({0.0 - 3 * dp, 1.0, 0.0 - 3 * dp}, {1.0 + 3 * dp, 1.0 + 3 * dp, 1.0 + 3 * dp});

    Box<3, double> xminplate({
                                 0.0 - 3 * dp,
                                 0.0,
                                 0.0 - 3 * dp,
                             },
                             {0.0, 1.0, 1.0 + 3 * dp});
    Box<3, double> xmaxplate({
                                 1.0,
                                 0.0,
                                 0.0 - 3 * dp,
                             },
                             {1.0 + 3 * dp, 1.0, 1.0 + 3 * dp});

    Box<3, double> zminplate({0.0 - 3 * dp, 0.0, 0.0 - 3 * dp}, {1.0 + 3 * dp, 1.0, 0.0});
    Box<3, double> zmaxplate({0.0 - 3 * dp, 0.0, 1.0}, {1.0 + 3 * dp, 1.0, 1.0 + 3 * dp});

    openfpm::vector<Box<3, double>> obstacle_and_bound_box;

    // obstacle_and_bound_box.add(zminplate);
    // obstacle_and_bound_box.add(zmaxplate);
    // obstacle_and_bound_box.add(xminplate);
    // obstacle_and_bound_box.add(xmaxplate);
    obstacle_and_bound_box.add(yminplate);
    obstacle_and_bound_box.add(ymaxplate);

    // obstacle_and_bound_box.add(obstacle1);
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

    vd.map();
    // Now that we fill the vector with particles
    ModelCustom md;
    vd.addComputationCosts(md);
    vd.getDecomposition().decompose();
    vd.map();
    vd.ghost_get<type, rho, Pressure, velocity>();

    vd.write_frame("output/Geometry", 0);

    auto NN = vd.getCellList(2 * H);
    // Evolve
    size_t write = 1;
    size_t it = 0;
    size_t it_reb = 0;
    double t = 0.0;

    double dt = DtMin;

    while (t <= t_end)
    {
        Vcluster<> &v_cl = create_vcluster();
        timer it_time;
        it_reb++;
        if (it_reb == 200)
        {
            vd.map();
            ModelCustom md;
            vd.addComputationCosts(md);
            vd.getDecomposition().decompose();
            if (v_cl.getProcessUnitID() == 0)
                std::cout << "REBALANCED " << std::endl;
        }
        vd.map();
        // Calculate pressure from the density

        double max_visc = 0.0;
        // TODO  ghost_get performance
        //       check which quantities are necesary
        vd.ghost_get<type, rho, rho_prev, vicous_force, drho, Pressure, velocity, force_p, v_pre, x_pre>();
        // Calc forces
        size_t pressureIteration = 0.0;
        size_t pressureIteration_min = 3.0;
        double rho_e_max = 0.00;
        double compressibility = 0.01;

        calc_forces_and_drho(vd, NN, max_visc);

        while (pressureIteration < pressureIteration_min || rho_e_max > compressibility)
        {
            rho_e_max = 0.0;
            // LOGFunction("predict_v_and_x");
            predict_v_and_x(vd, NN, dt);
            // LOGFunction("EqState_incompressible");
            // TODO  Domain sync
            //       check which quantities are necesary
            // vd.ghost_get<type, rho, rho_prev,vicous_force,drho,Pressure, velocity,force_p,v_pre,x_pre>();
            // vd.map();
            EqState_incompressible(vd, NN, max_visc, rho_e_max, dt);
            // LOGFunction("extrapolate_Boundaries");
            extrapolate_Boundaries(vd, NN, max_visc);
            // LOGFunction("calc_PressureForces");
            calc_PressureForces(vd, NN, max_visc);

            v_cl.max(rho_e_max);
            v_cl.execute();

            if ((rho_e_max > 0.01) && (v_cl.getProcessUnitID() == 0))
            {
                LOGdouble("pressureIteration : ", pressureIteration);
                LOGdouble("error max : ", rho_e_max);
            }
            ++pressureIteration;
        }

        // Get the maximum viscosity term across processors
        v_cl.max(max_visc);
        v_cl.execute();
        // Calculate delta t integration
        dt = calc_deltaT(vd, max_visc);
        // VerletStep or euler step
        /*
        it++;
        if (it < 40)
            peng_int(vd, dt);
        else
        {
            peng_int(vd, dt);
            it = 0;
        }
        */
        peng_int(vd, dt);
        it++;

        t += dt;

        if (it % 1 == 0)
        {

            // sensor_pressure calculation require ghost and update cell-list
            vd.map();
            vd.ghost_get<type, rho, Pressure, velocity>();
            vd.updateCellList(NN);
            // calculate the pressure at the sensor points

            vd.write_frame("output/Geometry", write);
            it_reb = 0;
            write++;

            if (v_cl.getProcessUnitID() == 0)
            {
                std::cout << "TIME: " << t << "  write " << it_time.getwct() << "   " << v_cl.getProcessUnitID() << "   " << cnt << "   Max visc: " << max_visc << std::endl;
            }
        }
        else
        {
            if (v_cl.getProcessUnitID() == 0)
            {
                std::cout << "TIME: " << t << "  " << it_time.getwct() << "   " << v_cl.getProcessUnitID() << "   " << cnt << "    Max visc: " << max_visc << std::endl;
            }
        }
    }
    openfpm_finalize();
}