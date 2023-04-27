// #define SE_CLASS1
// #define STOP_ON_ERROR
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include <limits>

#include "main.h"

int main(int argc, char *argv[])
{
    // initialize the library
    openfpm_init(&argc, &argv);
    // It contain for each time-step the value detected by the probes
    openfpm::vector<openfpm::vector<double>> press_t;
    openfpm::vector<Point<3, double>> probes;
    probes.add({0.8779, 0.3, 0.02});
    probes.add({0.754, 0.31, 0.02});
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
    // of Dual-SPH that use a different crpressureIterationia to draw particles
    Box<3, double> fluid_box({dp / 2.0, dp / 2.0, dp / 2.0}, {0.4 + dp / 2.0, 0.67 - dp / 2.0, 0.3 + dp / 2.0});
    // return an pressureIterationator to the fluid particles to add to vd
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);
    // here we fill some of the constants needed by the simulation
    max_fluid_height = fluid_it.getBoxMargins().getHigh(2);
    h_swl = fluid_it.getBoxMargins().getHigh(2) - fluid_it.getBoxMargins().getLow(2);
    B = (coeff_sound) * (coeff_sound)*gravity * h_swl * rho_zero / gamma_;
    cbar = coeff_sound * sqrt(gravity * h_swl);
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
    Box<3, double> recipient1({0.0 - 3 * dp, 0.0 - 3 * dp, 0.0 - 3.0 * dp}, {1.6 + dp * 3.0, 0.67 + dp * 3.0, 0.4 - dp / 2});
    // Box<3, double> recipient1({0.0, 0.0, 0.0}, {1.6 + dp / 2.0, 0.67 + dp / 2.0, 0.4 + dp / 2.0});
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
    // auto obstacle_box = DrawParticles::DrawSkin(vd, sz, domain, obstacle2, obstacle1);
    auto obstacle_box = DrawParticles::DrawBox(vd, sz, domain, obstacle1);
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
    vd.map();
    // Now that we fill the vector with particles
    ModelCustom md;
    vd.addComputationCosts(md);
    vd.getDecomposition().decompose();
    vd.map();
    vd.ghost_get<type, rho, Pressure, velocity>();
    auto NN = vd.getCellList(2 * H);
    // Evolve
    size_t write = 0;
    size_t it = 0;
    size_t it_reb = 0;
    double t = 0.0;
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
        //EqState(vd);
        double max_visc = 0.0;
        vd.ghost_get<type, rho, Pressure, velocity>();
        // Calc forces
        calc_nonPforces(vd, NN, max_visc);
        v_cl.max(max_visc);
        v_cl.execute();
        // Calculate delta t integration
        double dt = calc_deltaT(vd, max_visc);
        // Get the maximum viscosity term across processors
        double error_min = 0.01;
        int min_itteration = 3;
        int pressureIteration = 0;
        double rel_rho_predicted_error_max;

        // whiler max_error > something do
        while ((rel_rho_predicted_error_max > error_min) || (pressureIteration < min_itteration))
        {
            rel_rho_predicted_error_max = 0;
 
            predictPositionAndVelocity(vd, NN, dt, max_visc, rel_rho_predicted_error_max);
            EqState_incompressible(vd, NN, dt, max_visc, rel_rho_predicted_error_max);

            LOGExit("EqState", v_cl.getProcessUnitID());
            LOGEnter("entering  pressure force", v_cl.getProcessUnitID());

            extrapolate_Boundaries(vd, NN, max_visc);
            calc_forces_pressure(vd, NN, max_visc);

            LOGExit("leaving  pressure force", v_cl.getProcessUnitID());
            if (v_cl.getProcessUnitID() == 0)
            {
            LOGDouble("predicted error = ", rel_rho_predicted_error_max);
            LOGDouble("pressure pressureIterationation  = ", pressureIteration);
            LOGFunction("Leaving Pressure pressureIterationation");
            }

            ++pressureIteration;

            // if (pressureIteration > 25)
            // {
            //     break;
            // }
        }
        /* CHENG integration*/
        // VerletStep or euler step

        it++;
        if (it < 40)
            verlet_int(vd, dt);
        else
        {
            euler_int(vd, dt);
            it = 0;
        }
        

       //cheng_int(vd, dt);

        t += dt;
        if (it_reb == 50)
        {
            // sensor_pressure calculation require ghost and update cell-list
            vd.map();
            vd.ghost_get<type, rho, Pressure, velocity>();
            vd.updateCellList(NN);
            // calculate the pressure at the sensor points
            sensor_pressure(vd, NN, press_t, probes);
            vd.write_frame("output/Geometry", write);
            it_reb = 0;
            ++write;

            if (v_cl.getProcessUnitID() == 0)
            {
                std::cout << "TIME: " << t << "  write " << it_time.getwct() << "   " << v_cl.getProcessUnitID() << "   " << cnt << "   Max visc: " << max_visc << std::endl;
            }
        }
    }
    openfpm_finalize();
}