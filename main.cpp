
#include "main.h"
#include "PCISPHFunctions.h"

int main(int argc, char *argv[])
{

    // std::ofstream outFile;
    outFile.open("logfile.txt", std::ios_base::app);

    // initialize the library
    openfpm_init(&argc, &argv);
    // setup domain
    particles vd = setUpDomain();

    // Now that we fill the vector with particles

    ModelCustom md;
    vd.addComputationCosts(md);
    vd.getDecomposition().decompose();
    vd.map();
    vd.ghost_get<type, rho, Pressure, Pressure_prev, velocity, pressure_acc, viscous_acc>();
    auto NN = vd.getCellList(2 * H);

    // Decompose domain, map paricles and ghet ghost properties
    vd.write_frame("output/Geometry", 0);

    // Inti
    size_t write = 1;
    size_t nr_timestep = 0;
    size_t it_reb = 0;
    double t = 0.0;

    std::cout << "enter time loop " << std::endl;
    while (t <= t_end)
    {

        Vcluster<> &v_cl = create_vcluster();
        timer it_time;
        it_reb++;
        if (it_reb == 200)
        {
            vd.map();
            it_reb = 0;
            ModelCustom md;
            vd.addComputationCosts(md);
            vd.getDecomposition().decompose();
            if (v_cl.getProcessUnitID() == 0)
                std::cout << "REBALANCED " << std::endl;
        }
        vd.map();

        // initializeParameter
        double max_visc = 0.0;
        double densityError = 0.0;
        int pcisphIt = 0;
        const double dt = 0.05 * H / 3.0;

        vd.ghost_get<type, rho, Pressure, Pressure_prev, velocity, pressure_acc, viscous_acc>();

        Init_Loop(vd, NN, max_visc);

        vd.ghost_get<type, rho, Pressure, Pressure_prev, velocity, pressure_acc, viscous_acc>();
        
        while (densityError > maxDensityVariation || pcisphIt < 3)
        {
            densityError = 0.0;

            position_and_velocity_prediction(vd, NN, dt);
            predictDensity(vd, NN, max_visc, dt, densityError);
            predictPressure(vd, NN, max_visc, dt, densityError);
            resetPosVelDen(vd, NN);
            extrapolateBoundaries(vd, NN, max_visc);
            calc_Pressure_forces(vd, NN, max_visc);

            // get max densityError and max_art_visc
            v_cl.max(max_visc);
            v_cl.max(densityError);
            v_cl.execute();

            if (v_cl.getProcessUnitID() == 0 && ((pcisphIt > 0) && (pcisphIt % 5 == 0)))
            {
                outFile << "error :  " << densityError << std::endl;
            }
            ++pcisphIt;
        }

        stroemer_verlet_int(vd, dt);
        t += dt;
        nr_timestep++;

        printAndLog(vd, NN, write, cnt, max_visc, pcisphIt, t, it_reb, nr_timestep);
    }
    openfpm_finalize();
}