// #define SE_CLASS1
// #define STOP_ON_ERROR
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include "inputs.h"
#include "vanillaFunctions.h"
#include "UtilityFunctions.h"
#include "pengBauingerFunctions.h"

int main(int argc, char *argv[])
{
    // initialize the library
    openfpm_init(&argc, &argv);
    particles vd = InitializeDomain_PoiseuilleFlow3d();

    vd.map();
    vd.ghost_get<type, rho, Pressure, velocity>();
    vd.write_frame("output/Geometry", 0);
    auto NN = vd.getCellList(2 * H);
    // Evolve
    size_t write = 1;
    size_t it = 0;
    size_t it_reb = 0;
    double t = 0.0;
    while (t <= t_end)
    {
        Vcluster<> &v_cl = create_vcluster();
        timer it_time;
        it_reb++;

        vd.map();
        // Calculate pressure from the density

        EqState(vd);
        double max_visc = 0.0;

        vd.ghost_get<type, rho, Pressure, velocity>();
        calc_forces_peng(vd, NN, max_visc);
        // Get the maximum viscosity term across processors
        v_cl.max(max_visc);
        v_cl.execute();
 
        // Calculate delta t integration
        double dt = calc_deltaT(vd, max_visc);
        // VerletStep or euler step
        it++;

        peng_int(vd, dt);

        t += dt;
        if (write < t * 20)
        {
            // sensor_pressure calculation require ghost and update cell - list
            vd.map();
            vd.ghost_get<type, rho, Pressure, velocity>();

            vd.write_frame("output/Geometry", write);
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
