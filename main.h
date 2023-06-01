// #define SE_CLASS1
// #define STOP_ON_ERROR
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include <limits>
// A constant to indicate boundary particles
#define BOUNDARY 0
// A constant to indicate fluid particles
#define FLUID 1
// initial spacing between particles dp in the formulas
double eps = std::numeric_limits<double>::epsilon();

#define LOGFunction(x) std::cout << x << '\n'
#define LOGEnter(x, y) std::cout << " Entering function " << x << " from processor " << y << "\n"
#define LOGExit(x, y) std::cout << " Exiting function " << x << " from processor " << y << "\n"
#define LOGDouble(x, y) std::cout << " " << x << " " << y << "\n"

const double dp = 0.0085;
// Maximum height of the fluid water
// is going to be calculated and filled later on
double h_swl = 0.0;
// c_s in the formulas (constant used to calculate the sound speed)
double B = 0.0;
const double coeff_sound = 20.0;
// gamma in the formulas
const double gamma_ = 7.0;
// sqrt(3.0*dp*dp) support of the kernel
const double H = 0.0147224318643;
// Eta in the formulas
const double Eta2 = 0.01 * H * H;
// alpha in the formula
const double visco = 0.1;
// alpha in the formula
const double dynamic_viscosity = 1e-3;
// cbar in the formula (calculated later)
double cbar = 0.0;
// Mass of the fluid particles
const double MassFluid = 0.000614125;
// Mass of the boundary particles
const double MassBound = 0.000614125;
// End simulation time
#ifdef TEST_RUN
const double t_end = 0.001;
#else
const double t_end = 1.5;
#endif
// Gravity acceleration
const double gravity = 9.81;
// Reference densitu 1000Kg/m^3
const double rho_zero = 1000.0;

// Constant used to define time integration
const double CFLnumber = 0.2;
// Minimum T
const double DtMin =  0.05* H /1.9 ;
// Minimum Rho allowed
const double RhoMin = 700.0;
// Maximum Rho allowed
const double RhoMax = 1300.0;
// Filled in initialization
double max_fluid_height = 0.0;
// Properties
// FLUID or BOUNDARY
const size_t type = 0;
// Density
const int rho = 1;
// Density at step n-1
const int rho_prev = 2;
// Pressure
const int Pressure = 3;
// Delta rho calculated in the force calculation
const int drho = 4;
// calculated force
const int vicous_force = 5;
// velocity
const int velocity = 6;
// velocity at previous step
const int velocity_prev = 7;
// Type of the vector containing particles
const int x_pre = 8;
// Type of the vector containing particles
const int v_pre = 9;
// Type of the vector containing particles
const int rho_pred = 10;
// Type of the vector containing particles
const int force_p = 11;
// Type of the vector containing particles
const double bodyforce[3] = {0.0, 0.0, -9.81};

#define LOG(x) std::cout << x << '\n'

#define LOGdouble(x, y) std::cout << x << "  " << y << '\n'

typedef vector_dist<3, double, aggregate<size_t, double, double, double, double, double[3], double[3], double[3], double[3], double[3], double, double[3]>> particles;
//                                       |        |        |       |         |        |       |           |           |
//                                       |        |        |       |         |        |       |           |           |
//                                     type     density   density Pressure delta   vicous_force   velocity   velocityy   force
//                                                        at n-1           density                    at n - 1    pressure
struct ModelCustom
{
    template <typename Decomposition, typename vector>
    inline void addComputation(Decomposition &dec,
                               vector &vd,
                               size_t v,
                               size_t p)
    {
        if (vd.template getProp<type>(p) == FLUID)
            dec.addComputationCost(v, 4);
        else
            dec.addComputationCost(v, 3);
    }
    template <typename Decomposition>
    inline void applyModel(Decomposition &dec, size_t v)
    {
        dec.setSubSubDomainComputationCost(v, dec.getSubSubDomainComputationCost(v) * dec.getSubSubDomainComputationCost(v));
    }
    double distributionTol()
    {
        return 1.01;
    }
};
inline void EqState(particles &vd)
{
    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto a = it.get();
        double rho_a = vd.template getProp<rho>(a);
        double rho_frac = rho_a / rho_zero;
        vd.template getProp<Pressure>(a) = B * (rho_frac * rho_frac * rho_frac * rho_frac * rho_frac * rho_frac * rho_frac - 1.0);
        ++it;
    }
}
const double a2 = 1.0 / M_PI / H / H / H;
inline double Wab(double r)
{
    r /= H;
    if (r < 1.0)
        return (1.0 - 3.0 / 2.0 * r * r + 3.0 / 4.0 * r * r * r) * a2;
    else if (r < 2.0)
        return (1.0 / 4.0 * (2.0 - r) * (2.0 - r) * (2.0 - r)) * a2;
    else
        return 0.0;
}
const double c1 = -3.0 / M_PI / H / H / H / H;
const double d1 = 9.0 / 4.0 / M_PI / H / H / H / H;
const double c2 = -3.0 / 4.0 / M_PI / H / H / H / H;
const double a2_4 = 0.25 * a2;
// Filled later
double W_dap = 0.0;
inline void DWab(Point<3, double> &dx, Point<3, double> &DW, double r, bool print)
{
    const double qq = r / H;
    double qq2 = qq * qq;
    double fac1 = (c1 * qq + d1 * qq2) / r;
    double b1 = (qq < 1.0) ? 1.0f : 0.0f;
    double wqq = (2.0 - qq);
    double fac2 = c2 * wqq * wqq / r;
    double b2 = (qq >= 1.0 && qq < 2.0) ? 1.0f : 0.0f;
    double factor = (b1 * fac1 + b2 * fac2);
    DW.get(0) = factor * dx.get(0);
    DW.get(1) = factor * dx.get(1);
    DW.get(2) = factor * dx.get(2);
}
template <typename CellList>
inline void calc_forces_and_drho(particles &vd, CellList &NN, double &max_visc)
{
    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);
    while (part.isNext())
    {
        // ... a
        auto a = part.get();

        Point<3, double> xa = vd.getPos(a);
        Point<3, double> va = vd.getProp<velocity>(a);

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);
        double Va = (massa / rhoa);

        vd.template getProp<vicous_force>(a)[0] = 0.0;
        vd.template getProp<vicous_force>(a)[1] = 0.0;
        vd.template getProp<vicous_force>(a)[2] = 0.0;

        vd.template getProp<force_p>(a)[0] = 0.0;
        vd.template getProp<force_p>(a)[1] = 0.0;
        vd.template getProp<force_p>(a)[2] = 0.0;
        vd.template getProp<Pressure>(a) = 0.0;

        double v_force_x = 0.0;
        double v_force_y = 0.0;
        double v_force_z = 0.0;

        double viscosityaveraging = (2 * dynamic_viscosity * dynamic_viscosity) / (dynamic_viscosity + dynamic_viscosity);

        if (vd.getProp<type>(a) == FLUID)
        {
            auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
            while (Np.isNext() == true)
            {
                auto b = Np.get();
                Point<3, double> xb = vd.getPos(b);
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                };

                Point<3, double> vb = vd.getProp<velocity>(b);
                Point<3, double> dr = xa - xb;

                double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
                double rhob = vd.getProp<rho>(b);
                double Vb = (massa / rhoa);
                double r2 = norm2(dr);

                if (r2 < 4.0 * H * H)
                {
                    Point<3, double> v_rel = va - vb;
                    Point<3, double> DW;
                    double r = sqrt(r2);
                    Point<3, double> e_ab = dr / r;

                    Vb = (massb / rhob);

                    DWab(dr, DW, r, false);

                    double factor = (viscosityaveraging) * (Vb * Vb + Va * Va) * (DW.get(0) * e_ab.get(0) + DW.get(1) * e_ab.get(1) + DW.get(2) * e_ab.get(2));

                    v_force_x += factor * v_rel.get(0) / r;
                    v_force_y += factor * v_rel.get(1) / r;
                    v_force_z += factor * v_rel.get(2) / r;
                }
                ++Np;
            }
            vd.getProp<vicous_force>(a)[0] = v_force_x / massa * 0;
            vd.getProp<vicous_force>(a)[1] = v_force_y / massa * 0;
            vd.getProp<vicous_force>(a)[2] = v_force_z / massa * 0;
        }
        ++part;
    }
}
template <typename CellList>
inline void EqState_incompressible(particles &vd, CellList &NN, double &max_visc, double &rho_e_max, double &dt)
{
    auto part = vd.getDomainIterator();

    vd.updateCellList(NN);
    double density_pred_error;
    double rho_e;

    while (part.isNext())
    {
        // ... a
        auto a = part.get();

        Point<3, double> xa = vd.getPos(a);
        Point<3, double> va = vd.getProp<velocity>(a);
        Point<3, double> term_1_vec = {0, 0, 0};

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);

        if (vd.getProp<type>(a) == FLUID)
        {
            double pressureKoefficient = 0.0;
            double term_1_sca = 0.0;
            double term_2_sca = 0.0;
            double kernel = 0.0;
            double density_pred = 0.0;
            vd.getProp<drho>(a) = 0.0;

            auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
            while (Np.isNext() == true)
            {
                auto b = Np.get();
                Point<3, double> xb = vd.getPos(b);
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                };

                Point<3, double> vb = vd.getProp<velocity>(b);
                Point<3, double> dr = xa - xb;

                double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
                double rhob = vd.getProp<rho>(b);
                double r2 = norm2(dr);

                if (r2 < 4.0 * H * H)
                {
                    Point<3, double> v_rel = va - vb;
                    Point<3, double> DW;
                    double r = sqrt(r2);

                    DWab(dr, DW, r, false);
                    vd.getProp<drho>(a) += massb * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
                    term_1_vec += DW;
                    //TODO hier nur einmall masse und pcisph funktioniert
                    term_2_sca += (DW.get(0) * DW.get(0) + DW.get(1) * DW.get(1) + DW.get(2) * DW.get(2));
                }
                ++Np;
            };

            term_1_sca = term_1_vec.get(0) * term_1_vec.get(0) + term_1_vec.get(1) * term_1_vec.get(1) + term_1_vec.get(2) * term_1_vec.get(2);
            double beta = term_1_sca + term_2_sca;
            pressureKoefficient = (rho_zero * rho_zero) / (massa*dt * dt * beta);
            //TODO add comparisson using vd.getProp<rho_pred>(a) and assigning to real rho
            //  this would meand that the real densities are used for the calculations
            
            density_pred = vd.getProp<rho>(a) + vd.getProp<drho>(a) * dt;
            // vd.getProp<rho>(a) = vd.getProp<rho>(a) + vd.getProp<drho>(a) * dt;
            density_pred_error = density_pred - rho_zero;
            rho_e = std::abs((density_pred_error) / rho_zero);
            rho_e_max = std::max(rho_e_max, rho_e);

            double candidatePressure = vd.getProp<Pressure>(a)  + pressureKoefficient * density_pred_error  ;

            //TODO No Pressure Clipping 
            //  this would meand that the real densities are used for the calculations

            vd.getProp<Pressure>(a) =  candidatePressure    ;
        }

        ++part;
    }
}

template <typename CellList>
inline void calc_PressureForces(particles &vd, CellList &NN, double &max_visc)
{
    auto part = vd.getDomainIterator();
    vd.updateCellList(NN);
    while (part.isNext())
    {
        // ... a
        auto a = part.get();

        Point<3, double> xa = vd.getPos(a);
        Point<3, double> va = vd.getProp<velocity>(a);

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);
        double Va = (massa / rhoa);
        double Pa = vd.getProp<Pressure>(a);

        double p_force_x = 0.0;
        double p_force_y = 0.0;
        double p_force_z = 0.0;

        if (vd.getProp<type>(a) == FLUID)
        {
            auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
            while (Np.isNext() == true)
            {
                auto b = Np.get();
                Point<3, double> xb = vd.getPos(b);
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                };

                Point<3, double> vb = vd.getProp<velocity>(b);
                Point<3, double> dr = xa - xb;

                double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
                double Pb = vd.getProp<Pressure>(b);
                double rhob = vd.getProp<rho>(b);
                double Vb = (massa / rhoa);
                double r2 = norm2(dr);

                if (r2 < 4.0 * H * H)
                {
                    Point<3, double> v_rel = va - vb;
                    Point<3, double> DW;
                    double r = sqrt(r2);

                    Vb = (massb / rhob);
                    DWab(dr, DW, r, false);
                    double factor = (-1) * (Vb * Vb + Va * Va) * (rhob * vd.getProp<Pressure>(a) + rhoa * vd.getProp<Pressure>(b)) / (rhob + rhob);

                    p_force_x += factor * DW.get(0);
                    p_force_y += factor * DW.get(1);
                    p_force_z += factor * DW.get(2);
                }
                ++Np;
            }
            vd.getProp<force_p>(a)[0] = p_force_x / massa;
            vd.getProp<force_p>(a)[1] = p_force_y / massa;
            vd.getProp<force_p>(a)[2] = p_force_z / massa;
        }
        ++part;
    }
}
template <typename CellList>
inline void extrapolate_Boundaries(particles &vd, CellList &NN, double &max_visc)
{
    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        auto a = part.get();
        Point<3, double> xa = vd.getPos(a);
        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);
        double Va = (massa / rhoa);
        double Pa = vd.getProp<Pressure>(a);
        Point<3, double> va = vd.getProp<velocity>(a);
        if (vd.getProp<type>(a) != FLUID)
        {
            auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
            double kernel = 0;
            double rho_wab = 0;
            double p_wab = 0;
            double g_wab = 0;
            // For each neighborhood particle
            while (Np.isNext() == true)
            {
                auto b = Np.get();
                Point<3, double> xb = vd.getPos(b);
                double pressure_b = vd.getProp<Pressure>(b);

                if (a.getKey() == b || vd.getProp<type>(b) != FLUID)
                {
                    ++Np;
                    continue;
                };
                double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
                double rhob = vd.getProp<rho>(b);
                Point<3, double> dr = xa - xb;
                double r2 = norm2(dr);

                if (r2 < 4.0 * H * H)
                {
                    double wab = Wab(sqrt(r2));
                    kernel += wab;
                    rho_wab += rhob * wab;
                    p_wab += pressure_b * wab;
                    g_wab += (wab * rhob * bodyforce[0] * dr.get(0) + wab * rhob * bodyforce[1] * dr.get(1) + wab * rhob * bodyforce[2] * dr.get(2));
                }
                ++Np;
            }
            vd.getProp<rho>(a) = (rho_wab / kernel > rho_zero) ? rho_wab / kernel : rho_zero;
            vd.getProp<Pressure>(a) = ((p_wab + g_wab) / kernel > 0) ? (p_wab + g_wab) / kernel : 0.0;
        }
        ++part;
    }
}
void max_acceleration_and_velocity(particles &vd, double &max_acc, double &max_vel)
{
    // Calculate the maximum acceleration
    auto part = vd.getDomainIterator();
    while (part.isNext())
    {

        auto a = part.get();
        Point<3, double> acc(bodyforce);
        double acc2 = norm2(acc);
        Point<3, double> vel(vd.getProp<velocity>(a));
        double vel2 = norm2(vel);
        if (vel2 >= max_vel)
            max_vel = vel2;
        if (acc2 >= max_acc)
            max_acc = acc2;
        ++part;
    }
    max_acc = sqrt(max_acc);
    max_vel = sqrt(max_vel);
    Vcluster<> &v_cl = create_vcluster();
    v_cl.max(max_acc);
    v_cl.max(max_vel);
    v_cl.execute();
}
double calc_deltaT(particles &vd, double ViscDtMax)
{
    double Maxvel   = 1.9               ;
    double dt       = (0.05* H /1.9) ;
    
    // double Maxacc = 0.0;
    // double Maxvel = 0.0;
    // max_acceleration_and_velocity(vd, Maxacc, Maxvel);
    // //-dt1 depends on force per unit mass.
    // const double dt_f = (Maxacc) ? sqrt(H / Maxacc) : std::numeric_limits<int>::max();
    // //-dt2 combines the Courant and the viscous time-step controls.
    // const double dt_cv = H / (std::max(cbar, Maxvel * 10.) + H * ViscDtMax);
    // //-dt new value of time step.
    // double dt = double(CFLnumber) * std::min(dt_f, dt_cv);
    // if (dt < double(DtMin))
    //     dt = double(DtMin);
    // return dt;
    // starting with constant timestep

    return dt;
}
openfpm::vector<size_t> to_remove;
size_t cnt = 0;
void verlet_int(particles &vd, double dt)
{
    // list of the particle to remove
    to_remove.clear();
    // particle iterator
    auto part = vd.getDomainIterator();
    double dt205 = dt * dt * 0.5;
    double dt2 = dt * 2.0;
    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();
        // if the particle is boundary
        if (vd.template getProp<type>(a) == BOUNDARY)
        {

            // Update only the density
            vd.template getProp<velocity>(a)[0] = 0.0;
            vd.template getProp<velocity>(a)[1] = 0.0;
            vd.template getProp<velocity>(a)[2] = 0.0;
            // double rhonew = vd.template getProp<rho_prev>(a) + dt2*vd.template getProp<drho>(a);
            // vd.template getProp<rho>(a) = (rhonew < rho_zero)?rho_zero:rhonew;
            vd.template getProp<rho_prev>(a) = vd.template getProp<rho>(a);
            ;
            ++part;
            continue;
        }
        //-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
        double dx = vd.template getProp<velocity>(a)[0] * dt + (vd.template getProp<vicous_force>(a)[0] + vd.template getProp<force_p>(a)[0] + bodyforce[0]) * dt205;
        double dy = vd.template getProp<velocity>(a)[1] * dt + (vd.template getProp<vicous_force>(a)[1] + vd.template getProp<force_p>(a)[1] + bodyforce[1]) * dt205;
        double dz = vd.template getProp<velocity>(a)[2] * dt + (vd.template getProp<vicous_force>(a)[2] + vd.template getProp<force_p>(a)[2] + bodyforce[2]) * dt205;
        vd.getPos(a)[0] += dx;
        vd.getPos(a)[1] += dy;
        vd.getPos(a)[2] += dz;
        double velX = vd.template getProp<velocity>(a)[0];
        double velY = vd.template getProp<velocity>(a)[1];
        double velZ = vd.template getProp<velocity>(a)[2];
        double rhop = vd.template getProp<rho>(a);
        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + (vd.template getProp<vicous_force>(a)[0] + vd.template getProp<force_p>(a)[0] + bodyforce[0]) * dt2;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + (vd.template getProp<vicous_force>(a)[1] + vd.template getProp<force_p>(a)[1] + bodyforce[1]) * dt2;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + (vd.template getProp<vicous_force>(a)[2] + vd.template getProp<force_p>(a)[2] + bodyforce[2]) * dt2;

        double rho_candidate = vd.template getProp<rho_prev>(a) + dt2 * vd.template getProp<drho>(a);
        vd.template getProp<rho>(a) = (rho_candidate > rho_zero) ? rho_candidate : rho_zero;
        // vd.template getProp<rho>(a) = vd.template getProp<rho_prev>(a) + dt2 * vd.template getProp<drho>(a);

        // Check if the particle go out of range in space and in density
        if (vd.getPos(a)[0] < 0.000263878 || vd.getPos(a)[1] < 0.000263878 || vd.getPos(a)[2] < 0.000263878 ||
            vd.getPos(a)[0] > 0.000263878 + 1.59947 || vd.getPos(a)[1] > 0.000263878 + 0.672972 || vd.getPos(a)[2] > 0.000263878 + 0.903944 ||
            vd.template getProp<rho>(a) < RhoMin || vd.template getProp<rho>(a) > RhoMax)
        {
            to_remove.add(a.getKey());
        }
        vd.template getProp<velocity_prev>(a)[0] = velX;
        vd.template getProp<velocity_prev>(a)[1] = velY;
        vd.template getProp<velocity_prev>(a)[2] = velZ;
        vd.template getProp<rho_prev>(a) = rhop;
        ++part;
    }
    // remove the particles
    vd.remove(to_remove, 0);
    // increment the iteration counter
    cnt++;
}
void euler_int(particles &vd, double dt)
{
    // list of the particle to remove
    to_remove.clear();
    // particle iterator
    auto part = vd.getDomainIterator();
    double dt205 = dt * dt * 0.5;
    double dt2 = dt * 2.0;
    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();
        // if the particle is boundary
        if (vd.template getProp<type>(a) == BOUNDARY)
        {
            double rhop = vd.template getProp<rho>(a);
            // Update only the density

            vd.template getProp<velocity>(a)[1] = 0.0;
            vd.template getProp<velocity>(a)[2] = 0.0;
            ++part;
            continue;
        }
        //-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
        double dx = vd.template getProp<velocity>(a)[0] * dt + (vd.template getProp<vicous_force>(a)[0] + vd.template getProp<force_p>(a)[0] + bodyforce[0]) * dt205;
        double dy = vd.template getProp<velocity>(a)[1] * dt + (vd.template getProp<vicous_force>(a)[1] + vd.template getProp<force_p>(a)[1] + bodyforce[1]) * dt205;
        double dz = vd.template getProp<velocity>(a)[2] * dt + (vd.template getProp<vicous_force>(a)[2] + vd.template getProp<force_p>(a)[2] + bodyforce[2]) * dt205;

        vd.getPos(a)[0] = vd.template getProp<x_pre>(a)[0] + dx;
        vd.getPos(a)[1] = vd.template getProp<x_pre>(a)[1] + dy;
        vd.getPos(a)[2] = vd.template getProp<x_pre>(a)[2] + dz;

        vd.template getProp<x_pre>(a)[0] = vd.getPos(a)[0];
        vd.template getProp<x_pre>(a)[1] = vd.getPos(a)[1];
        vd.template getProp<x_pre>(a)[2] = vd.getPos(a)[2];

        double velX = vd.template getProp<velocity>(a)[0];
        double velY = vd.template getProp<velocity>(a)[1];
        double velZ = vd.template getProp<velocity>(a)[2];
        double rhop = vd.template getProp<rho>(a);

        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity>(a)[0] + (vd.template getProp<vicous_force>(a)[0] + vd.template getProp<force_p>(a)[0] + bodyforce[0]) * dt;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity>(a)[1] + (vd.template getProp<vicous_force>(a)[1] + vd.template getProp<force_p>(a)[1] + bodyforce[1]) * dt;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity>(a)[2] + (vd.template getProp<vicous_force>(a)[2] + vd.template getProp<force_p>(a)[2] + bodyforce[2]) * dt;

        double rho_candidate = vd.template getProp<rho_prev>(a) + dt2 * vd.template getProp<drho>(a);
        vd.template getProp<rho>(a) = (rho_candidate > rho_zero) ? rho_candidate : rho_zero;
        // vd.template getProp<rho>(a) = vd.template getProp<rho_prev>(a) + dt2 * vd.template getProp<drho>(a);

        // Check if the particle go out of range in space and in density
        if (vd.getPos(a)[0] < 0.000263878 || vd.getPos(a)[1] < 0.000263878 || vd.getPos(a)[2] < 0.000263878 ||
            vd.getPos(a)[0] > 0.000263878 + 1.59947 || vd.getPos(a)[1] > 0.000263878 + 0.672972 || vd.getPos(a)[2] > 0.000263878 + 0.903944 ||
            vd.template getProp<rho>(a) < RhoMin || vd.template getProp<rho>(a) > RhoMax)
        {
            std::cout << "removing stuff" << '\n';
            to_remove.add(a.getKey());
        }

        vd.template getProp<velocity_prev>(a)[0] = velX;
        vd.template getProp<velocity_prev>(a)[1] = velY;
        vd.template getProp<velocity_prev>(a)[2] = velZ;
        vd.template getProp<rho_prev>(a) = rhop;
        ++part;
    }
    // remove the particles
    vd.remove(to_remove, 0);
    // increment the iteration counter
    cnt++;
}
void peng_int(particles &vd, double dt)
{
    // list of the particle to remove
    to_remove.clear();
    // particle iterator
    auto part = vd.getDomainIterator();
    double dt05 = dt * 0.5;
    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();
        // if the particle is boundary
        if (vd.template getProp<type>(a) == BOUNDARY)
        {
            // double rhop = vd.template getProp<rho>(a);
            //  Update only the density
            vd.template getProp<velocity>(a)[0] = 0.0;
            vd.template getProp<velocity>(a)[1] = 0.0;
            vd.template getProp<velocity>(a)[2] = 0.0;
            ++part;
            continue;
        }

        //-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.

        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + (vd.template getProp<vicous_force>(a)[0] + vd.template getProp<force_p>(a)[0] + bodyforce[0]) * dt;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + (vd.template getProp<vicous_force>(a)[1] + vd.template getProp<force_p>(a)[1] + bodyforce[1]) * dt;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + (vd.template getProp<vicous_force>(a)[2] + vd.template getProp<force_p>(a)[2] + bodyforce[2]) * dt;

        vd.getPos(a)[0] = vd.template getProp<x_pre>(a)[0] + (vd.template getProp<velocity_prev>(a)[0] + vd.template getProp<velocity>(a)[0]) * dt05;
        vd.getPos(a)[1] = vd.template getProp<x_pre>(a)[1] + (vd.template getProp<velocity_prev>(a)[1] + vd.template getProp<velocity>(a)[1]) * dt05;
        vd.getPos(a)[2] = vd.template getProp<x_pre>(a)[2] + (vd.template getProp<velocity_prev>(a)[2] + vd.template getProp<velocity>(a)[2]) * dt05;

        //  TODO reached timestep is going to be wrong because of drho,
        //  proper handling needs to be assured
        double rho_candidate = vd.template getProp<rho_prev>(a) + dt * vd.template getProp<drho>(a);
        vd.template getProp<rho>(a) = (rho_candidate > rho_zero) ? rho_candidate : rho_zero;
        // vd.template getProp<rho>(a) = vd.template getProp<rho_prev>(a) + dt2 * vd.template getProp<drho>(a);

        // Check if the particle go out of range in space and in density
        if (vd.getPos(a)[0] < 0.000263878 || vd.getPos(a)[1] < 0.000263878 || vd.getPos(a)[2] < 0.000263878 ||
            vd.getPos(a)[0] > 0.000263878 + 1.59947 || vd.getPos(a)[1] > 0.000263878 + 0.672972 || vd.getPos(a)[2] > 0.000263878 + 0.903944 ||
            vd.template getProp<rho>(a) < RhoMin || vd.template getProp<rho>(a) > RhoMax)
        {
            std::cout << "removing stuff" << '\n';
            to_remove.add(a.getKey());
        }

        vd.template getProp<velocity_prev>(a)[0] = vd.template getProp<velocity>(a)[0];
        vd.template getProp<velocity_prev>(a)[1] = vd.template getProp<velocity>(a)[1];
        vd.template getProp<velocity_prev>(a)[2] = vd.template getProp<velocity>(a)[2];

        vd.template getProp<x_pre>(a)[0] = vd.getPos(a)[0];
        vd.template getProp<x_pre>(a)[1] = vd.getPos(a)[1];
        vd.template getProp<x_pre>(a)[2] = vd.getPos(a)[2];

        vd.template getProp<rho_prev>(a) = vd.template getProp<rho>(a);

        ++part;
    }
    // remove the particles
    vd.remove(to_remove, 0);
    // increment the iteration counter
    cnt++;
}

template <typename Vector, typename CellList>
inline void sensor_pressure(Vector &vd,
                            CellList &NN,
                            openfpm::vector<openfpm::vector<double>> &press_t,
                            openfpm::vector<Point<3, double>> &probes)
{
    Vcluster<> &v_cl = create_vcluster();
    press_t.add();
    for (size_t i = 0; i < probes.size(); i++)
    {
        float press_tmp = 0.0f;
        float tot_ker = 0.0;
        // if the probe is inside the processor domain
        if (vd.getDecomposition().isLocal(probes.get(i)) == true)
        {
            // Get the position of the probe i
            Point<3, double> xp = probes.get(i);
            // get the iterator over the neighbohood particles of the probes position
            auto itg = NN.template getNNIterator<NO_CHECK>(NN.getCell(probes.get(i)));
            while (itg.isNext())
            {
                auto q = itg.get();
                // Only the fluid particles are importants
                if (vd.template getProp<type>(q) != FLUID)
                {
                    ++itg;
                    continue;
                }
                // Get the position of the neighborhood particle q
                Point<3, double> xq = vd.getPos(q);
                // Calculate the contribution of the particle to the pressure
                // of the probe
                double r = sqrt(norm2(xp - xq));
                double ker = Wab(r) * (MassFluid / rho_zero);
                // Also keep track of the calculation of the summed
                // kernel
                tot_ker += ker;
                // Add the total pressure contribution
                press_tmp += vd.template getProp<Pressure>(q) * ker;
                // next neighborhood particle
                ++itg;
            }
            // We calculate the pressure normalizing the
            // sum over all kernels
            if (tot_ker == 0.0)
                press_tmp = 0.0;
            else
                press_tmp = 1.0 / tot_ker * press_tmp;
        }
        // This is not necessary in principle, but if you
        // want to make all processor aware of the history of the calculated
        // pressure we have to execute this
        v_cl.sum(press_tmp);
        v_cl.execute();
        // We add the calculated pressure into the history
        press_t.last().add(press_tmp);
    }
}
template <typename CellList>
void predict_v_and_x(particles &vd, CellList &NN, double dt)
{
    // list of the particle to remove
    to_remove.clear();
    // particle iterator
    auto part = vd.getDomainIterator();
    double dt05 = dt * 0.5;
    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();
        // if the particle is boundary
        if (vd.template getProp<type>(a) == BOUNDARY)
        {
            ++part;
            continue;
        }

        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + (vd.template getProp<vicous_force>(a)[0] + vd.template getProp<force_p>(a)[0] + bodyforce[0]) * dt;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + (vd.template getProp<vicous_force>(a)[1] + vd.template getProp<force_p>(a)[1] + bodyforce[1]) * dt;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + (vd.template getProp<vicous_force>(a)[2] + vd.template getProp<force_p>(a)[2] + bodyforce[2]) * dt;

        // vd.getPos(a)[0] = vd.template getProp<x_pre>(a)[0] + (vd.template getProp<velocity_prev>(a)[0] + vd.template getProp<velocity>(a)[0]) * dt05;
        // vd.getPos(a)[1] = vd.template getProp<x_pre>(a)[1] + (vd.template getProp<velocity_prev>(a)[1] + vd.template getProp<velocity>(a)[1]) * dt05;
        // vd.getPos(a)[2] = vd.template getProp<x_pre>(a)[2] + (vd.template getProp<velocity_prev>(a)[2] + vd.template getProp<velocity>(a)[2]) * dt05;

        //TODO prediction angepasst 
        
        vd.getPos(a)[0] = vd.template getProp<x_pre>(a)[0] + (vd.template getProp<velocity>(a)[0]) * dt;
        vd.getPos(a)[1] = vd.template getProp<x_pre>(a)[1] + (vd.template getProp<velocity>(a)[1]) * dt;
        vd.getPos(a)[2] = vd.template getProp<x_pre>(a)[2] + (vd.template getProp<velocity>(a)[2]) * dt;

        ++part;
    }
}