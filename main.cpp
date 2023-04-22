//#define SE_CLASS1
//#define STOP_ON_ERROR
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
// A constant to indicate boundary particles
#define BOUNDARY 0
// A constant to indicate fluid particles
#define FLUID 1


#include "main.h"
#include <chrono>
#include <vector>



//std::vector<float> time_per_timestep	;
//std::chrono::high_resolution_clock::time_point start, stop;

int main(int argc, char* argv[])
{

    openfpm_init(&argc,&argv);


    openfpm::vector<openfpm::vector<double>> press_t;
    openfpm::vector<Point<3,double>> probes;

    probes.add({0.8779,0.3,0.02});
    probes.add({0.754,0.31,0.02});

    Box<3,double> domain({-0.05,-0.05,-0.05},{1.7010,0.7065,0.5025});
    size_t sz[3] = {207,90,66};

    // Fill W_dap
    W_dap = 1.0/Wab(H/1.5);

    // Here we define the boundary conditions of our problem
    size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

    // extended boundary around the domain, and the processor domain
    Ghost<3,double> g(2*H);
    
    particles vd(0,domain,bc,g,DEC_GRAN(512));
    Box<3,double> fluid_box({dp/2.0,dp/2.0,dp/2.0},{0.4+dp/2.0,0.67-dp/2.0,0.3+dp/2.0});

    // return an iterator to the fluid particles to add to vd
    auto fluid_it = DrawParticles::DrawBox(vd,sz,domain,fluid_box);

    // here we fill some of the constants needed by the simulation
    max_fluid_height = fluid_it.getBoxMargins().getHigh(2);
    h_swl = fluid_it.getBoxMargins().getHigh(2) - fluid_it.getBoxMargins().getLow(2);
    B = (coeff_sound)*(coeff_sound)*gravity*h_swl*rho_zero / gamma_;
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
        vd.template getLastProp<Pressure>() = rho_zero * gravity *  (max_fluid_height - fluid_it.get().get(2));

        vd.template getLastProp<rho>() = pow(vd.template getLastProp<Pressure>() / B + 1, 1.0/gamma_) * rho_zero;
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
    Box<3,double> recipient1({0.0,0.0,0.0},{1.6+dp/2.0,0.67+dp/2.0,0.4+dp/2.0});
    Box<3,double> recipient2({dp,dp,dp},{1.6-dp/2.0,0.67-dp/2.0,0.4+dp/2.0});

    Box<3,double> obstacle1({0.9,0.24-dp/2.0,0.0},{1.02+dp/2.0,0.36,0.45+dp/2.0});
    Box<3,double> obstacle2({0.9+dp,0.24+dp/2.0,0.0},{1.02-dp/2.0,0.36-dp,0.45-dp/2.0});
    Box<3,double> obstacle3({0.9+dp,0.24,0.0},{1.02,0.36,0.45});

    openfpm::vector<Box<3,double>> holes;
    holes.add(recipient2);
    holes.add(obstacle1);
    auto bound_box = DrawParticles::DrawSkin(vd,sz,domain,holes,recipient1);

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
    auto obstacle_box = DrawParticles::DrawSkin(vd,sz,domain,obstacle2,obstacle1);

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

	//! \cond [load balancing] \endcond

    vd.ghost_get<type,rho,Pressure,velocity>();

    auto NN = vd.getCellList(2*H);

    size_t write = 0;
    size_t it = 0;
    size_t it_reb = 0;
    double t = 0.0;
	
	double dt =	DtMin;

    while (t <= t_end)
    {
		//start = std::chrono::high_resolution_clock::now();

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

        double max_visc = 0.0;

        vd.ghost_get<type,rho,Pressure,velocity>();

		// Calc all forces but pressure
		//calc_forces_withoutPressure(vd,NN,max_visc);

        LOGEnter("Calc_forces",v_cl.getProcessUnitID())     ;
        calc_forces(vd,NN,max_visc) ;
        LOGExit("Calc_forces",v_cl.getProcessUnitID())      ;
        /*
        LOGEnter("calc_Density",v_cl.getProcessUnitID())     ;
        */
        //  calc_Density(vd, NN)                            ;
        /*
        LOGExit("calc_Density",v_cl.getProcessUnitID())      ;
        */
		EqState_incompressible(vd,NN,dt,max_visc);
        LOGEnter("EqState",v_cl.getProcessUnitID())     ;
        EqState(vd)             ;
        LOGExit("EqState",v_cl.getProcessUnitID())      ;

        // Get the maximum viscosity term across processors
        v_cl.max(max_visc);
        v_cl.execute();

        // Calculate delta t integration
        double dt = calc_deltaT(vd,max_visc);

        // VerletStep or euler step
        it++;
        if (it < 40)
            verlet_int(vd,dt);
            //  verlet_int_no_densityUpdate(vd,dt);
        else
        {
            euler_int(vd,dt);
            //euler_intno_densityUpdate(vd,dt);
            it = 0;
        }

        t += dt;
		//if (write < t*100)
        if (true)
        {
            // sensor_pressure calculation require ghost and update cell-list
            vd.map();
            vd.ghost_get<type,rho,Pressure,velocity>();
            vd.updateCellList(NN);
            // calculate the pressure at the sensor points
            //sensor_pressure(vd,NN,press_t,probes);
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
		// stop = std::chrono::high_resolution_clock::now()		;
		// float microseconds = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
		// time_per_timestep.push_back(microseconds)	;
		// std::cout<< "microseconds per loop" << '\n'			;
		// std::cout<< microseconds << '\n'
    }

    openfpm_finalize();
}
 
