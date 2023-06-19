
#include "main.h"

int main(int argc, char *argv[])
{
    // initialize the library
    openfpm_init(&argc, &argv);
    // It contain for each time-step the value detected by the probes
    openfpm::vector<openfpm::vector<double>> press_t;
    openfpm::vector<Point<3, double>> probes;

    double Xmin = 0;
    double Ymin = 0 - 3 * dp;
    double Zmin = -3 * dp;

    double Xmax = 1.0;
    double Ymax = 1.0 + 3 * dp;
    double Zmax = +3 * dp;

    double L_x = Xmax - Xmin;
    double L_y = Ymax - Ymin;
    double L_z = Zmax - Zmin;

    size_t Nr_x = ceil(L_x / dp);
    size_t Nr_y = ceil(L_y / dp);
    size_t Nr_z = ceil(L_z / dp);

    std::cout << "Nr of particles in x, y, z : " << Nr_x << " " << Nr_y << " " << Nr_z << " "
              << "\n";

    // Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
    Box<3, double> domain({Xmin, Ymin, Zmin}, {Xmax, Ymax, Zmax});

    size_t sz[3] = {Nr_x, Nr_y, Nr_z};
    // Fill W_dap
    W_dap = 1.0 / Wab(H / 1.5);
    // Here we define the boundary conditions of our problem
    size_t bc[3] = {PERIODIC, NON_PERIODIC, PERIODIC};
    // extended boundary around the domain, and the processor domain
    Ghost<3, double> g(2 * H);

    particles vd(0, domain, bc, g, DEC_GRAN(512));

    Box<3, double> fluid_box({Xmin + 0.1 * dp, 0.0, Zmin + 0.1 * dp}, {Xmax - 0.1 * dp, 1.0, Zmax - 0.1 * dp});
    // return an iterator to the fluid particles to add to vd
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);
    // here we fill some of the constants needed by the simulation
    max_fluid_height = fluid_it.getBoxMargins().getHigh(2);
    h_swl = fluid_it.getBoxMargins().getHigh(2) - fluid_it.getBoxMargins().getLow(2);
    B = (coeff_sound) * (coeff_sound)*gravity * h_swl * rho_zero / gamma_;
    cbar = coeff_sound * sqrt(gravity * h_swl);
    int NrOfFluidParticles = 0;
    while (fluid_it.isNext())
    {
        vd.add();
        vd.getLastPos()[0] = fluid_it.get().get(0);
        vd.getLastPos()[1] = fluid_it.get().get(1);
        vd.getLastPos()[2] = fluid_it.get().get(2);

        vd.template getLastProp<type>() = FLUID;
        vd.template getLastProp<Pressure>() = rho_zero * gravity * (max_fluid_height - fluid_it.get().get(2)) + backgroundPressure;
        vd.template getLastProp<rho>() = pow(vd.template getLastProp<Pressure>() / B + 1, 1.0 / gamma_) * rho_zero;
        vd.template getLastProp<rho_prev>() = vd.template getLastProp<rho>();
        vd.template getLastProp<velocity>()[0] = 0.0;
        vd.template getLastProp<velocity>()[1] = 0.0;
        vd.template getLastProp<velocity>()[2] = 0.0;
        vd.template getLastProp<velocity_prev>()[0] = 0.0;
        vd.template getLastProp<velocity_prev>()[1] = 0.0;
        vd.template getLastProp<velocity_prev>()[2] = 0.0;

        // vd.template getLastProp<viscousFoces>()[0] = 0.0;
        // vd.template getLastProp<viscousFoces>()[1] = 0.0;
        // vd.template getLastProp<viscousFoces>()[2] = 0.0;

        ++fluid_it;
        ++NrOfFluidParticles;
    }

    Box<3, double> ymin({Xmin + 0.1 * dp, 0.0 - 3 * dp, Zmin + 0.1 * dp}, {Xmax - 0.1 * dp, 0.00, Zmax - 0.1 * dp});
    Box<3, double> ymax({Xmin + 0.1 * dp, 1.0, Zmin + 0.1 * dp}, {Xmax - 0.1 * dp, 1.0 + dp * 3.0, Zmax - 0.1 * dp});
    // obstacle box
    // Box<3, double> obstacle1({0.9, 0.24 - dp / 2.0, 0.0}, {1.02 + dp / 2.0, 0.36, 0.5});

    openfpm::vector<Box<3, double>> obstacle_and_bound_box;

    obstacle_and_bound_box.add(ymin);
    obstacle_and_bound_box.add(ymax);
    // obstacle_and_bound_box.add(obstacle1);

    std::cout << "draw box" << '\n';
    for (size_t i = 0; i < 2; ++i)
    {
        // std::cout << "for box number " << i << '\n';
        Box<3, double> box = obstacle_and_bound_box.get(i);
        auto toFill = DrawParticles::DrawBox(vd, sz, domain, box);

        while (toFill.isNext())
        {
            vd.add();
            vd.getLastPos()[0] = toFill.get().get(0);
            vd.getLastPos()[1] = toFill.get().get(1);
            vd.getLastPos()[2] = toFill.get().get(2);

            vd.template getLastProp<type>() = BOUNDARY;
            vd.template getLastProp<Pressure>() = backgroundPressure;
            vd.template getLastProp<rho>() = rho_zero;
            vd.template getLastProp<rho_prev>() = rho_zero;
            vd.template getLastProp<velocity>()[0] = 0.0;
            vd.template getLastProp<velocity>()[1] = 0.0;
            vd.template getLastProp<velocity>()[2] = 0.0;

            vd.template getLastProp<velocity_prev>()[0] = 0.0;
            vd.template getLastProp<velocity_prev>()[1] = 0.0;
            vd.template getLastProp<velocity_prev>()[2] = 0.0;

            // vd.template getLastProp<viscousFoces>()[0] = 0.0;
            // vd.template getLastProp<viscousFoces>()[1] = 0.0;
            // vd.template getLastProp<viscousFoces>()[2] = 0.0;

            ++toFill;
        }
    }

    // std::cout << "exit drawing box" << '\n';
    vd.map();
    vd.write_frame("output/Geometry", 0);
    // Now that we fill the vector with particles
    ModelCustom md;
    vd.addComputationCosts(md);
    vd.getDecomposition().decompose();
    vd.map();
    vd.ghost_get<type, rho, Pressure, velocity>();
    auto NN = vd.getCellList(2 * H);
    // Evolve
    size_t write_coounter = 1;
    size_t write = 1;
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
        // std::cout << "EqState " << std::endl;
        EqState(vd);
        double max_visc = 0.0;
        vd.ghost_get<type, rho, Pressure, velocity>();
        // Calc forces
        // std::cout << "enter calc_forces " << std::endl;
        calc_forces(vd, NN, max_visc);
        // std::cout << "exit calc_forces " << std::endl;
        // Get the maximum viscosity term across processors
        v_cl.max(max_visc);
        v_cl.execute();
        // Calculate delta t integration
        double dt = calc_deltaT(vd, max_visc);
        // VerletStep or euler step
        it++;
        write_coounter++;
        if (it < 40)
            verlet_int(vd, dt);
        else
        {
            euler_int(vd, dt);
            it = 0;
        }
        t += dt;
        if (write_coounter % 100 == 0)
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