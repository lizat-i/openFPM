
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
    // Decompose domain, map paricles and ghet ghost properties
    vd.map();
    vd.write_frame("output/Geometry", 0);
    vd.getDecomposition().decompose();
    vd.map();
    vd.ghost_get<>();
    auto NN = vd.getCellList(2 * H);

    // Inti
    size_t write = 1;
    size_t it = 0;
    double t = 0.0;

    while (t <= t_end)
    {
        Vcluster<> &v_cl = create_vcluster();

        vd.map();
        vd.template ghost_get<>();

        // initializeParameter
        double max_visc = 0.0;
        double densityError = 0.0;
        int iterCount = 0;
        const double dt = 0.05 * H / 3.0;

        Init_Loop(vd, NN, max_visc);

        while (densityError > maxDensityVariation || iterCount < 3)
        {
            densityError = 0.0;

            position_and_velocity_prediction(vd, dt);
            predictDensity(vd, NN, max_visc, dt, densityError);
            predictPressure(vd, NN, max_visc, dt, densityError);
            extrapolateBoundaries(vd, NN, max_visc);
            calc_Pressure_forces(vd, NN, max_visc);

            // get max densityError and max_art_visc
            v_cl.max(max_visc);
            v_cl.max(densityError);
            v_cl.execute();
            ++iterCount;
        }

        stroemer_verlet_int(vd, dt);
        t += dt;
        it++;

        printAndLog(vd, NN, write, cnt, max_visc, iterCount, t);
    }
    openfpm_finalize();
}