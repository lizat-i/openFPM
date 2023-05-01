// #define SE_CLASS1
// #define STOP_ON_ERROR
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include <limits>

#include "main.h"

int main(int argc, char *argv[])
{
    std::cout<<"started programm"<<'\n';


    openfpm_init(&argc, &argv);
    openfpm::vector<openfpm::vector<double>> press_t;
    openfpm::vector<Point<3, double>> probes;
    probes.add({0.8779, 0.3, 0.02});
    probes.add({0.754, 0.31, 0.02});
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
    cbar = coeff_sound * sqrt(gravity * h_swl);


    while (fluid_it.isNext())
    {
        vd.add();
        vd.getLastPos()[0] = fluid_it.get().get(0);
        vd.getLastPos()[1] = fluid_it.get().get(1);
        vd.getLastPos()[2] = fluid_it.get().get(2);
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
    }
    
    // Build domain with 6 boxes.
    // first groundplate
    Box<3, double> groundplate({0.0 - 3 * dp, 0.0 - 3 * dp, 0.0 - 3.0 * dp}, {1.6 + dp * 3.0, 0.67 + dp * 3.0, 0.0});
    // xmin xmax plates
    Box<3, double> xmin({0.0 - 3 * dp, 0.0 - 3 * dp, 0.0}, {0.0, 0.67 + dp * 3.0, 0.45});
    Box<3, double> xmax({1.6, 0.0 - 3 * dp, 0.0}, {1.6 + dp * 3.0, 0.67 + dp * 3.0, 0.45});
    // ymin ymax plates
    Box<3, double> ymin({0.0 - 3 * dp, 0.0 - 3 * dp, 0.0}, {1.6 + dp * 3.0, 0.00, 0.4});
    Box<3, double> ymax({0.0 - 3 * dp, 0.67, 0.0}, {1.6 + dp * 3.0, 0.67 + dp * 3.0, 0.4});
    // obstacle box
    Box<3, double> obstacle1({0.9, 0.24 - dp / 2.0, 0.0}, {1.02 + dp / 2.0, 0.36, 0.45});

    openfpm::vector<Box<3, double>> obstacle_and_bound_box;

    obstacle_and_bound_box.add(groundplate);
    obstacle_and_bound_box.add(xmin);
    obstacle_and_bound_box.add(xmax);
    obstacle_and_bound_box.add(ymin);
    obstacle_and_bound_box.add(ymax);
    obstacle_and_bound_box.add(obstacle1);
    std::cout<<"draw box"<<'\n';
    for (size_t i = 0; i < 6; ++i)
    {
        std::cout<<"for box number "<< i <<'\n';
        Box<3, double> box = obstacle_and_bound_box.get(i);
        auto toFill = DrawParticles::DrawBox(vd, sz, domain,box);

        while (toFill.isNext())
        {
            vd.add();
            vd.getLastPos()[0] = toFill.get().get(0);
            vd.getLastPos()[1] = toFill.get().get(1);
            vd.getLastPos()[2] = toFill.get().get(2);
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
        EqState(vd);
        double max_visc = 0.0;
        vd.ghost_get<type, rho, Pressure, velocity>();
        // Calc forces
        calc_forces(vd, NN, max_visc);
        // Get the maximum viscosity term across processors
        v_cl.max(max_visc);
        v_cl.execute();
        // Calculate delta t integration
        double dt = calc_deltaT(vd, max_visc);
        // VerletStep or euler step
        it++;
        if (it < 40)
            verlet_int(vd, dt);
        else
        {
            euler_int(vd, dt);
            it = 0;
        }
        t += dt;
 
        // sensor_pressure calculation require ghost and update cell-list
        vd.map();
        vd.ghost_get<type, rho, Pressure, velocity>();
        vd.updateCellList(NN);
        // calculate the pressure at the sensor points
        sensor_pressure(vd, NN, press_t, probes);
        vd.write_frame("output/Geometry", write);
        it_reb = 0;
        //++write;

        if (v_cl.getProcessUnitID() == 0)
        {
            std::cout << "TIME: " << t << "  write " << it_time.getwct() << "   " << v_cl.getProcessUnitID() << "   " << cnt << "   Max visc: " << max_visc << std::endl;
        }
    }
    openfpm_finalize();
}