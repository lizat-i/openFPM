//#define SE_CLASS1
//#define STOP_ON_ERROR
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include "inputs.h"
#include "vanillaFunctions.h"
#include "UtilityFunctions.h"

int main(int argc, char* argv[])
{
    // initialize the library
    openfpm_init(&argc,&argv);

    particles vd = InitializeDomain();

    vd.map();
    // Now that we fill the vector with particles
    ModelCustom md;
    vd.addComputationCosts(md);
    vd.getDecomposition().decompose();
    vd.map();
    vd.ghost_get<type,rho,Pressure,velocity>();
    auto NN = vd.getCellList(2*H);
    // Evolve
    size_t write = 0;
    size_t it = 0;
    size_t it_reb = 0;
    double t = 0.0;
    while (t <= t_end)
    {
        Vcluster<> & v_cl = create_vcluster();
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
        // Calculate pressure from the density
        EqState(vd);
        double max_visc = 0.0;
        vd.ghost_get<type,rho,Pressure,velocity>();
        // Calc forces
        calc_forces(vd,NN,max_visc);
        // Get the maximum viscosity term across processors
        v_cl.max(max_visc);
        v_cl.execute();
        // Calculate delta t integration
        double dt = calc_deltaT(vd,max_visc);
        // VerletStep or euler step
        it++;
        if (it < 40)
            verlet_int(vd,dt);
        else
        {
            euler_int(vd,dt);
            it = 0;
        }
        t += dt;
        if (write < t*100)
        {
            // sensor_pressure calculation require ghost and update cell-list
            vd.map();
            vd.ghost_get<type,rho,Pressure,velocity>();
 
            vd.write_frame("output/Geometry",write);
            write++;
            if (v_cl.getProcessUnitID() == 0)
            {std::cout << "TIME: " << t << "  write " << it_time.getwct() << "   " << v_cl.getProcessUnitID() << "   " << cnt << "   Max visc: " << max_visc << std::endl;}
        }
        else
        {
            if (v_cl.getProcessUnitID() == 0)
            {std::cout << "TIME: " << t << "  " << it_time.getwct() << "   " << v_cl.getProcessUnitID() << "   " << cnt << "    Max visc: " << max_visc << std::endl;}
        }
    }
    openfpm_finalize();
}
 