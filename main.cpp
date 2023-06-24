
#include "main.h"

int main(int argc, char *argv[])
{

    // std::ofstream outFile;
    outFile.open("logfile.txt", std::ios_base::app);

    // initialize the library
    openfpm_init(&argc, &argv);
    // It contain for each time-step the value detected by the probes
    openfpm::vector<openfpm::vector<double>> press_t;
    openfpm::vector<Point<3, double>> probes;

    Xmin = 0 - 3 * dp;
    Ymin = 0 - 3 * dp;
    Zmin = 0 - 3 * dp;

    Xmax = 1.6 + 3 * dp;
    Ymax = 0.67 + 3 * dp;
    Zmax = 0.45;

    L_x = Xmax - Xmin;
    L_y = Ymax - Ymin;
    L_z = Zmax - Zmin;

    size_t Nr_x = ceil(L_x / dp);
    size_t Nr_y = ceil(L_y / dp);
    size_t Nr_z = ceil(L_z / dp);

    std::cout << "Nr of particles in x, y, z : " << Nr_x << " " << Nr_y << " " << Nr_z << " "
              << "\n";
    Box<3, double> domain({Xmin, Ymin, Zmin}, {Xmax, Ymax, Zmax});
    size_t sz[3] = {Nr_x, Nr_y, Nr_z};

    W_dap = 1.0 / Wab(H / 1.5);
    size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
    Ghost<3, double> g(2 * H);
    particles vd(0, domain, bc, g, DEC_GRAN(512));
    Box<3, double> fluid_box({0.0, 0.0, 0.0}, {0.4, 0.67, 0.3});
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

    // here we fill some of the constants needed by the simulation
    max_fluid_height = fluid_it.getBoxMargins().getHigh(2);
    h_swl = fluid_it.getBoxMargins().getHigh(2) - fluid_it.getBoxMargins().getLow(2);
    B = (coeff_sound) * (coeff_sound)*gravity * h_swl * rho_zero / gamma_;
    cbar = coeff_sound * sqrt(gravity * h_swl);
    cbar = 3.0 *10 ; 
    std::cout<<cbar<< std::endl;

    int NrOfFluidParticles = 0;
    while (fluid_it.isNext())
    {
        vd.add();
        vd.getLastPos()[0] = fluid_it.get().get(0);
        vd.getLastPos()[1] = fluid_it.get().get(1);
        vd.getLastPos()[2] = fluid_it.get().get(2);
        vd.getLastProp<x_pre>()[0] = fluid_it.get().get(0);
        vd.getLastProp<x_pre>()[1] = fluid_it.get().get(1);
        vd.getLastProp<x_pre>()[2] = fluid_it.get().get(2);
        vd.template getLastProp<type>() = FLUID;
        vd.template getLastProp<Pressure>() = 0     ;   //rho_zero * gravity * (max_fluid_height - fluid_it.get().get(2));
        vd.template getLastProp<rho>() = rho_zero   ;   //pow(vd.template getLastProp<Pressure>() / B + 1, 1.0 / gamma_) * rho_zero;
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

    // Build domain with 6 boxes.
    // first groundplate
    Box<3, double> groundplate({0.0 - 3 * dp, 0.0 - 3 * dp, 0.0 - 3.0 * dp}, {1.6 + dp * 3.0, 0.67 + dp * 3.0, 0.0});

    Box<3, double> xmin({0.0 - 3 * dp, 0.0 - 3 * dp, 0.0}, {0.0, 0.67 + dp * 3.0, 0.40});
    Box<3, double> xmax({1.6 - dp / 2, 0.0 - 3 * dp, 0.0}, {1.6 + dp * 3.0, 0.67 + dp * 3.0, 0.4});

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

    std::cout << "draw box" << '\n';
    for (size_t i = 0; i < 5; ++i)
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

            vd.getLastProp<x_pre>()[0] = toFill.get().get(0);;
            vd.getLastProp<x_pre>()[1] = toFill.get().get(1);;
            vd.getLastProp<x_pre>()[2] = toFill.get().get(2);;

            vd.template getLastProp<type>() = BOUNDARY;
            vd.template getLastProp<Pressure>() = 0;
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
        // Domain Housekeepng
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

        // Enter relevant loop
        double max_visc = 0.0;
        Init_Loop(vd, NN, max_visc);

        double densityError = 0.0;
        int iterCount = 0;
        const double dt = 0.05 * H / 3.0;

        while (densityError > maxDensityVariation || iterCount < 3)
        {
            densityError = 0.0;

            position_and_velocity_prediction(vd, dt);
            vd.map();
            vd.ghost_get<type, rho, Pressure, velocity>();
            vd.updateCellList(NN);

            predictDensityAndUpdate(vd, NN, max_visc, dt, densityError);
            // outFile << "enterExBoun " << iterCount << std::endl;
            extrapolateBoundaries(vd, NN, max_visc);
            // outFile << "entercalcPre " << iterCount << std::endl;
            calc_Pressure_forces(vd, NN, max_visc);

            v_cl.max(max_visc);
            v_cl.max(densityError);
            v_cl.execute();

            // v_cl.max(densityError);
            // v_cl.execute();

            ++iterCount;
        }
        if (v_cl.getProcessUnitID() == 0)
        {
            outFile << "it " << iterCount << std::endl;
        }

        it++;
        write_coounter++;
        stroemer_verlet_int(vd, dt);
        t += dt;
        if (write_coounter % 1 == 0)
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