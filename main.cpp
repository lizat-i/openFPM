#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include "main.h"

operationParams parameters = {
    /*BOUNDARY=*/0,
    /*FLUID=*/1,
    /*dp=*/0.0085,
    /*h_swl=*/0.0,
    /*coeff_sound=*/20.0,
    /*gamma_=*/7.0,
    /*H=*/0.0147224318643,
    /*Eta2=*/0.01 * 0.0147224318643 * 0.0147224318643,
    /*visco=*/0.1,
    /*cbar=*/0.0,
    /*MassFluid=*/0.000614125,
    /*MassBound=*/0.000614125,
    /*t_end=*/1.5,
    /*gravity=*/9.81,
    /*rho_zero=*/1000.0,
    /*B=*/0.0,
    /*W_dap = */ 0.0,
    /*CFLnumber=*/0.2,
    /*DtMin=*/0.00001,
    /*RhoMin=*/700.0,
    /*RhoMax=*/1300.0,
    /*max_fluid_height=*/0.0,
    /*a2=*/1.0 / M_PI / 0.0147224318643 / 0.0147224318643 / 0.0147224318643,
    /*c1=*/-3.0 / M_PI / 0.0147224318643 / 0.0147224318643 / 0.0147224318643 / 0.0147224318643,
    /*d1=*/9.0 / 4.0 / M_PI / 0.0147224318643 / 0.0147224318643 / 0.0147224318643 / 0.0147224318643,
    /*c2=*/-3.0 / 4.0 / M_PI / 0.0147224318643 / 0.0147224318643 / 0.0147224318643 / 0.0147224318643,
    /*a2_4=*/0.25 * (1.0 / M_PI / 0.0147224318643 / 0.0147224318643 / 0.0147224318643),
    /*cnt = */0,
    /*max_visc*/
};

int main(int argc, char *argv[])
{

  // initialize the library Domain and set Probes
  openfpm_init(&argc, &argv);
  openfpm::vector<openfpm::vector<double>> press_t;
  openfpm::vector<Point<3, double>> probes;
  probes.add({0.8779, 0.3, 0.02});
  probes.add({0.754, 0.31, 0.02});
  Box<3, double> domain({-0.05, -0.05, -0.05}, {1.7010, 0.7065, 0.5025});
  size_t sz[3] = {207, 90, 66};
  parameters.W_dap = 1.0 / Wab(parameters.H / 1.5,parameters);
  size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
  Ghost<3, double> g(2 * parameters.a2);
  particles vd(0, domain, bc, g, DEC_GRAN(512));


  createBoxAndParseBox(vd, parameters,domain,g);
  vd.map();

  vd.write("myParticles", VTK_WRITER);
  openfpm_finalize();
 
}