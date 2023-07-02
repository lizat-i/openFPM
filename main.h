// #define SE_CLASS1
// #define STOP_ON_ERROR
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include <limits>
#include <iostream>
#include <fstream>
#include "constants.h"

typedef vector_dist<3, double, aggregate<size_t, double, double, double, double, double[3], double[3], double[3], double[3], double[3], double>> particles;
//                                        |      |        |          |            |            |         |            |           |
//                                        |      |        |          |            |            |         |            |           |
//                                      type   density   density    Pressure    delta       force      velocity    velocity    position
//                                                       at n-1                 density                             at n - 1    at n - 1

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
        if (vd.template getProp<type>(a) == BOUNDARY)
        {
            ++it;
            continue;
        }
        double rho_a = vd.template getProp<rho>(a);
        double rho_frac = rho_a / rho_zero;
        double candidatePressure = B * (rho_frac * rho_frac * rho_frac * rho_frac * rho_frac * rho_frac * rho_frac - 1.0);
        vd.template getProp<Pressure>(a) = (candidatePressure > 0.0) ? candidatePressure : 0.0;
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
    double fac1 = (c1 * qq + d1 * qq2) / (r + EPS);
    double b1 = (qq < 1.0) ? 1.0f : 0.0f;
    double wqq = (2.0 - qq);
    double fac2 = c2 * wqq * wqq / (r + EPS);
    double b2 = (qq >= 1.0 && qq < 2.0) ? 1.0f : 0.0f;
    double factor = (b1 * fac1 + b2 * fac2);
    DW.get(0) = factor * dx.get(0);
    DW.get(1) = factor * dx.get(1);
    DW.get(2) = factor * dx.get(2);
}
// Tensile correction
inline double Tensile(double r, double rhoa, double rhob, double prs1, double prs2)
{
    const double qq = r / H;
    //-Cubic Spline kernel
    double wab;
    if (r > H)
    {
        double wqq1 = 2.0f - qq;
        double wqq2 = wqq1 * wqq1;
        wab = a2_4 * (wqq2 * wqq1);
    }
    else
    {
        double wqq2 = qq * qq;
        double wqq3 = wqq2 * qq;
        wab = a2 * (1.0f - 1.5f * wqq2 + 0.75f * wqq3);
    }
    //-Tensile correction.
    double fab = wab * W_dap;
    fab *= fab;
    fab *= fab; // fab=fab^4
    const double tensilp1 = (prs1 / (rhoa * rhoa)) * (prs1 > 0 ? 0.01 : -0.2);
    const double tensilp2 = (prs2 / (rhob * rhob)) * (prs2 > 0 ? 0.01 : -0.2);
    return (fab * (tensilp1 + tensilp2));
}
inline double Pi(const Point<3, double> &dr, double rr2, Point<3, double> &dv, double rhoa, double rhob, double massb, double &visc)
{
    const double dot = dr.get(0) * dv.get(0) + dr.get(1) * dv.get(1) + dr.get(2) * dv.get(2);
    const double dot_rr2 = dot / (rr2 + Eta2);
    visc = std::max(dot_rr2, visc);
    if (dot < 0)
    {
        const float amubar = H * dot_rr2;
        const float robar = (rhoa + rhob) * 0.5f;
        const float pi_visc = (-visco * cbar * amubar / robar);
        return pi_visc;
    }
    else
        return 0.0;
}
template <typename CellList>
inline void calc_forces(particles &vd, CellList &NN, double &max_visc)
{
    auto part = vd.getDomainIterator();

    // Update the cell-list
    // std::cout << "exit updateCellList" << std::endl;
    vd.updateCellList(NN);
    // std::cout << "exit updateCellList" << std::endl;

    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();

        Point<3, double> xa = vd.getPos(a);
        Point<3, double> va = vd.getProp<velocity>(a);

        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;
        double rhoa = vd.getProp<rho>(a);
        double Pa = vd.getProp<Pressure>(a);

        // Reset the force counter (- gravity on zeta direction)
        vd.template getProp<pressure_acc>(a)[0] = 0.0;
        vd.template getProp<pressure_acc>(a)[1] = 0.0;
        vd.template getProp<pressure_acc>(a)[2] = 0.0;

        vd.template getProp<viscous_acc>(a)[0] = 0.0;
        vd.template getProp<viscous_acc>(a)[1] = 0.0;
        vd.template getProp<viscous_acc>(a)[2] = 0.0;

        vd.template getProp<drho>(a) = 0.0;

        // We threat FLUID particle differently from BOUNDARY PARTICLES ...
        if (vd.getProp<type>(a) == BOUNDARY)
        {
            auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
            double kernel = 0;
            double rho_wab = 0;
            double p_wab = 0;
            double g_wab = 0;
            double v_wab = 0;

            Point<3, double> g_wab_vectorial = {0, 0, 0};
            Point<3, double> v_wab_vectorial = {0, 0, 0};
            // For each neighborhood particle
            while (Np.isNext() == true)
            {
                // ... q
                auto b = Np.get();

                // Get the position xp of the particle
                Point<3, double> xb = vd.getPos(b);
                Point<3, double> vb = vd.getProp<velocity>(b);

                double Pb = vd.getProp<Pressure>(b);

                if (a.getKey() == b || vd.getProp<type>(b) != FLUID)
                {
                    ++Np;
                    continue;
                };

                // get the mass of the particle
                double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
                double rhob = vd.getProp<rho>(b);
                Point<3, double> dr = xa - xb;

                double r2 = norm2(dr);
                double r = sqrt(r2);

                // If the particles interact ...
                if (r < smoothingRadius)
                {

                    double wab = Wab(r);
                    kernel += wab;
                    rho_wab += rhob * wab;
                    p_wab += Pb * wab;
                    v_wab_vectorial += vb * wab;
                    g_wab_vectorial += wab * rhob * dr;
                    // Pressure and density calclation Boundary Particle
                }

                ++Np;
            }
            g_wab = g_wab_vectorial.get(0) * bodyforce[0] + g_wab_vectorial.get(1) * bodyforce[1] + g_wab_vectorial.get(2) * bodyforce[2];

            double cand_pressure = (p_wab + g_wab) / kernel;
            double cand_density = rho_wab / kernel;

            vd.getProp<rho>(a) = (kernel > 0.0) ? cand_density : rho_zero;
            vd.getProp<Pressure>(a) = (kernel > 0.0) ? cand_pressure : 0.0;
            // vd.getProp<Pressure>(a) = cand_pressure ;

            vd.template getProp<velocity>(a)[0] = (kernel > 0.0) ? -v_wab_vectorial[0] / kernel : 0;
            vd.template getProp<velocity>(a)[1] = (kernel > 0.0) ? -v_wab_vectorial[1] / kernel : 0;
            vd.template getProp<velocity>(a)[2] = (kernel > 0.0) ? -v_wab_vectorial[2] / kernel : 0;
        }
        else
        {

            auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
            double viscousFocesTermx = 0.0;
            double viscousFocesTermy = 0.0;
            double viscousFocesTermz = 0.0;
            double kernel = 0;

            // For each neighborhood particle
            while (Np.isNext() == true)
            {
                // ... q
                auto b = Np.get();
                Point<3, double> xb = vd.getPos(b);
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                };

                Point<3, double> vb = vd.getProp<velocity>(b);
                Point<3, double> dr = xa - xb;
                double Pb = vd.getProp<Pressure>(b);
                double rhob = vd.getProp<rho>(b);
                double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
                double r2 = norm2(dr);
                double r = sqrt(r2);

                // If the particles interact ...
                if (r < smoothingRadius)
                {

                    Point<3, double> v_rel = va - vb;
                    Point<3, double> DW;
                    Point<3, double> e_ab = dr * (1 / r);
                    DWab(dr, DW, r, false);
                    double wab = Wab(r);
                    double Va = massa / rhoa;
                    double Vb = massb / rhob;
                    kernel += wab;

                    double factor = -massb * ((vd.getProp<Pressure>(a) + vd.getProp<Pressure>(b)) / (rhoa * rhob)); // + Tensile(r, rhoa, rhob, Pa, Pb) );

                    vd.getProp<pressure_acc>(a)[0] += factor * DW.get(0);
                    vd.getProp<pressure_acc>(a)[1] += factor * DW.get(1);
                    vd.getProp<pressure_acc>(a)[2] += factor * DW.get(2);

                    vd.getProp<drho>(a) += massb * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));

                    double viscousForcesFactor = -massb * Pi(dr, r2, v_rel, rhoa, rhob, massb, max_visc);

                    viscousFocesTermx += viscousForcesFactor * DW.get(0);
                    viscousFocesTermy += viscousForcesFactor * DW.get(1);
                    viscousFocesTermz += viscousForcesFactor * DW.get(2);

                    vd.getProp<drho>(a) += massb * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
                }
                ++Np;
            }

            vd.getProp<viscous_acc>(a)[0] = (kernel > 0.0) ? viscousFocesTermx : 0.0;
            vd.getProp<viscous_acc>(a)[1] = (kernel > 0.0) ? viscousFocesTermy : 0.0;
            vd.getProp<viscous_acc>(a)[2] = (kernel > 0.0) ? viscousFocesTermz : 0.0;
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
    double Maxacc = 0.0;
    double Maxvel = 0.0;
    max_acceleration_and_velocity(vd, Maxacc, Maxvel);
    //-dt1 depends on force per unit mass.
    const double dt_f = (Maxacc) ? sqrt(H / Maxacc) : std::numeric_limits<int>::max();
    //-dt2 combines the Courant and the viscous time-step controls.
    const double dt_cv = H / (std::max(cbar, Maxvel * 10.) + H * ViscDtMax);
    //-dt new value of time step.
    double dt = double(delta_t_coeeficient) * std::min(dt_f, dt_cv);
    if (dt < double(DtMin))
        dt = double(DtMin);
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
            // Update rho
            double rhop = vd.template getProp<rho>(a);
            // Update only the density
            // vd.template getProp<velocity>(a)[0] = 0.0;
            // vd.template getProp<velocity>(a)[1] = 0.0;
            // vd.template getProp<velocity>(a)[2] = 0.0;
            // double rhonew = vd.template getProp<rho_prev>(a) + dt2*vd.template getProp<drho>(a);
            // vd.template getProp<rho>(a) = (rhonew < rho_zero)?rho_zero:rhonew;
            vd.template getProp<rho_prev>(a) = rhop;
            ++part;
            continue;
        }
        //-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
        double dx = vd.template getProp<velocity>(a)[0] * dt + (vd.template getProp<pressure_acc>(a)[0] + vd.template getProp<viscous_acc>(a)[0] + bodyforce[0]) * dt205;
        double dy = vd.template getProp<velocity>(a)[1] * dt + (vd.template getProp<pressure_acc>(a)[1] + vd.template getProp<viscous_acc>(a)[1] + bodyforce[1]) * dt205;
        double dz = vd.template getProp<velocity>(a)[2] * dt + (vd.template getProp<pressure_acc>(a)[2] + vd.template getProp<viscous_acc>(a)[2] + bodyforce[2]) * dt205;
        vd.getPos(a)[0] += dx;
        vd.getPos(a)[1] += dy;
        vd.getPos(a)[2] += dz;
        double velX = vd.template getProp<velocity>(a)[0];
        double velY = vd.template getProp<velocity>(a)[1];
        double velZ = vd.template getProp<velocity>(a)[2];
        double rhop = vd.template getProp<rho>(a);

        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + (vd.template getProp<pressure_acc>(a)[0] + vd.template getProp<viscous_acc>(a)[0] + bodyforce[0]) * dt2;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + (vd.template getProp<pressure_acc>(a)[1] + vd.template getProp<viscous_acc>(a)[1] + bodyforce[1]) * dt2;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + (vd.template getProp<pressure_acc>(a)[2] + vd.template getProp<viscous_acc>(a)[2] + bodyforce[2]) * dt2;

        double CandidateDensity = vd.template getProp<rho_prev>(a) + dt2 * vd.template getProp<drho>(a);
        vd.template getProp<rho>(a) = (CandidateDensity > rho_zero) ? CandidateDensity : rho_zero;

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
            // Update rho
            double rhop = vd.template getProp<rho>(a);
            // Update only the density
            // vd.template getProp<velocity>(a)[0] = 0.0;
            // vd.template getProp<velocity>(a)[1] = 0.0;
            // vd.template getProp<velocity>(a)[2] = 0.0;
            // double rhonew = vd.template getProp<rho>(a) + dt*vd.template getProp<drho>(a);
            // vd.template getProp<rho>(a) = (rhonew < rho_zero)?rho_zero:rhonew;
            ++part;
            continue;
        }
        //-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
        double dx = vd.template getProp<velocity>(a)[0] * dt + (vd.template getProp<pressure_acc>(a)[0] + vd.template getProp<viscous_acc>(a)[0] + bodyforce[0]) * dt205;
        double dy = vd.template getProp<velocity>(a)[1] * dt + (vd.template getProp<pressure_acc>(a)[1] + vd.template getProp<viscous_acc>(a)[1] + bodyforce[1]) * dt205;
        double dz = vd.template getProp<velocity>(a)[2] * dt + (vd.template getProp<pressure_acc>(a)[2] + vd.template getProp<viscous_acc>(a)[2] + bodyforce[2]) * dt205;
        vd.getPos(a)[0] += dx;
        vd.getPos(a)[1] += dy;
        vd.getPos(a)[2] += dz;

        double velX = vd.template getProp<velocity>(a)[0];
        double velY = vd.template getProp<velocity>(a)[1];
        double velZ = vd.template getProp<velocity>(a)[2];
        double rhop = vd.template getProp<rho>(a);
        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity>(a)[0] + (vd.template getProp<pressure_acc>(a)[0] + vd.template getProp<viscous_acc>(a)[0] + bodyforce[0]) * dt;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity>(a)[1] + (vd.template getProp<pressure_acc>(a)[1] + vd.template getProp<viscous_acc>(a)[1] + bodyforce[1]) * dt;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity>(a)[2] + (vd.template getProp<pressure_acc>(a)[2] + vd.template getProp<viscous_acc>(a)[2] + bodyforce[2]) * dt;

        double CandidateDensity = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);
        vd.template getProp<rho>(a) = (CandidateDensity > rho_zero) ? CandidateDensity : rho_zero;

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
inline void printAndLog(particles &vd, CellList &NN, size_t &write, size_t &cnt, double &max_visc, int &pcisphIt, double &t, size_t it_reb, size_t nr_timestep)
{
    Vcluster<> &v_cl = create_vcluster();
    auto part = vd.getDomainIterator();

    // Log Itercount Print other stuff
    if (v_cl.getProcessUnitID() == 0)
    {
        outFile << "pcisphIt :  " << pcisphIt << std::endl;
    }
    if (nr_timestep % 5 == 0)
    {
        // sensor_pressure calculation require ghost and update cell-list
        vd.map();
        vd.ghost_get<type, rho, Pressure, Pressure_prev, velocity, pressure_acc, viscous_acc>();
        vd.updateCellList(NN);

        vd.write_frame("output/Geometry", write);
        ++write;

        if (v_cl.getProcessUnitID() == 0)
        {
            std::cout << "TIME: " << t << " process unit " << v_cl.getProcessUnitID() << "   " << cnt << "   Max visc: " << max_visc << std::endl;
        }
    }
}

inline particles setUpDomain()
{

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
    cbar = 3.0 * 10;
    std::cout << cbar << std::endl;

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
        vd.template getLastProp<Pressure>() = 0; // rho_zero * gravity * (max_fluid_height - fluid_it.get().get(2));
        vd.template getLastProp<Pressure_prev>() = 0;
        vd.template getLastProp<rho>() = rho_zero; // pow(vd.template getLastProp<Pressure>() / B + 1, 1.0 / gamma_) * rho_zero;
        vd.template getLastProp<rho_prev>() = vd.template getLastProp<rho>();
        vd.template getLastProp<velocity>()[0] = 0.0;
        vd.template getLastProp<velocity>()[1] = 0.0;
        vd.template getLastProp<velocity>()[2] = 0.0;
        vd.template getLastProp<velocity_prev>()[0] = 0.0;
        vd.template getLastProp<velocity_prev>()[1] = 0.0;
        vd.template getLastProp<velocity_prev>()[2] = 0.0;
        vd.template getLastProp<rho_err>() = 0.0;

        // vd.template getLastProp<viscous_acc>()[0] = 0.0;
        // vd.template getLastProp<viscous_acc>()[1] = 0.0;
        // vd.template getLastProp<viscous_acc>()[2] = 0.0;

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

            vd.getLastProp<x_pre>()[0] = toFill.get().get(0);
            vd.getLastProp<x_pre>()[1] = toFill.get().get(1);
            vd.getLastProp<x_pre>()[2] = toFill.get().get(2);

            vd.template getLastProp<type>() = BOUNDARY;
            vd.template getLastProp<Pressure>() = 0;
            vd.template getLastProp<Pressure_prev>() = 0;
            vd.template getLastProp<rho>() = rho_zero;

            vd.template getLastProp<rho_prev>() = rho_zero;
            vd.template getLastProp<velocity>()[0] = 0.0;
            vd.template getLastProp<velocity>()[1] = 0.0;
            vd.template getLastProp<velocity>()[2] = 0.0;

            vd.template getLastProp<velocity_prev>()[0] = 0.0;
            vd.template getLastProp<velocity_prev>()[1] = 0.0;
            vd.template getLastProp<velocity_prev>()[2] = 0.0;
            vd.template getLastProp<rho_err>() = 0.0;
            // vd.template getLastProp<viscous_acc>()[0] = 0.0;
            // vd.template getLastProp<viscous_acc>()[1] = 0.0;
            // vd.template getLastProp<viscous_acc>()[2] = 0.0;

            ++toFill;
        }
    }

    return vd;
}
inline double dotProduct(Point<3, double> &a, Point<3, double> &b)
{

    return a.get(0) * b.get(0) + a.get(1) * b.get(1) + a.get(2) * b.get(2);
}