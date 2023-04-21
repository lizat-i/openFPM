//! \cond [inclusion] \endcond
#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"

#define LOGFunction(x) // std::cout << x << '\n'

// initial spacing between particles dp in the formulas
const double dp = 0.0085;
// Maximum height of the fluid water
// is going to be calculated and filled later on
double h_swl = 0.0;

// c_s in the formulas (constant used to calculate the sound speed)
const double coeff_sound = 20.0;

// gamma in the formulas
const double gamma_ = 7.0;

// sqrt(3.0*dp*dp) support of the kernel
const double H = 0.0147224318643;

// Eta in the formulas
const double Eta2 = 0.01 * H * H;

// alpha in the formula
const double visco = 0.1;

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
// const double t_end = 1.5;
const double t_end = 0.005;
#endif

// Gravity acceleration
const double gravity = 9.81;

// Reference densitu 1000Kg/m^3
const double rho_zero = 1000.0;

// Filled later require h_swl, it is b in the formulas
double B = 0.0;

// Constant used to define time integration
const double CFLnumber = 0.2;

// Minimum T
const double DtMin = 0.00001;

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
const int force = 5;

// velocity
const int velocity = 6;

// velocity at previous step
const int velocity_prev = 7;
// Type of the vector containing particles
typedef vector_dist<3, double, aggregate<size_t, double, double, double, double, double[3], double[3], double[3]>> particles;
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
    LOGFunction(" a2 ");
    LOGFunction(a2);
    LOGFunction("q");
    LOGFunction(r);

    LOGFunction("Output Wab \n");

	if (r < 1.0)
		return (1.0 - 3.0 / 2.0 * r * r + 3.0 / 4.0 * r * r * r) * a2;
	else if (r < 2.0)
		return (1.0 / 4.0 * (2.0 - r) * (2.0 - r ) * (2.0 - r)) * a2;
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
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();

        // Get the position xp of the particle
        Point<3, double> xa = vd.getPos(a);

        // Take the mass of the particle dependently if it is FLUID or BOUNDARY
        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;

        // Get the density of the of the particle a
        double rhoa = vd.getProp<rho>(a);

        // Get the pressure of the particle a
        double Pa = vd.getProp<Pressure>(a);

        // Get the Velocity of the particle a
        Point<3, double> va = vd.getProp<velocity>(a);

        // Reset the force counter (- gravity on zeta direction)
        vd.template getProp<force>(a)[0] = 0.0;
        vd.template getProp<force>(a)[1] = 0.0;
        vd.template getProp<force>(a)[2] = -gravity;
        vd.template getProp<drho>(a) = 0.0;

        // We threat FLUID particle differently from BOUNDARY PARTICLES ...
        if (vd.getProp<type>(a) != FLUID)
        {
            // If it is a boundary particle calculate the delta rho based on equation 2
            // This require to run across the neighborhoods particles of a
            auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

            // For each neighborhood particle
            while (Np.isNext() == true)
            {
                // ... q
                auto b = Np.get();

                // Get the position xp of the particle
                Point<3, double> xb = vd.getPos(b);

                // if (p == q) skip this particle
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                };

                // get the mass of the particle
                double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;

                // Get the velocity of the particle b
                Point<3, double> vb = vd.getProp<velocity>(b);

                // Get the pressure and density of particle b
                double Pb = vd.getProp<Pressure>(b);
                double rhob = vd.getProp<rho>(b);

                // Get the distance between p and q
                Point<3, double> dr = xa - xb;
                // take the norm of this vector
                double r2 = norm2(dr);

                // If the particles interact ...
                if (r2 < 4.0 * H * H)
                {
                    // ... calculate delta rho
                    double r = sqrt(r2);

                    Point<3, double> dv = va - vb;

                    Point<3, double> DW;
                    DWab(dr, DW, r, false);

                    const double dot = dr.get(0) * dv.get(0) + dr.get(1) * dv.get(1) + dr.get(2) * dv.get(2);
                    const double dot_rr2 = dot / (r2 + Eta2);
                    max_visc = std::max(dot_rr2, max_visc);

                    vd.getProp<drho>(a) += massb * (dv.get(0) * DW.get(0) + dv.get(1) * DW.get(1) + dv.get(2) * DW.get(2));
                }

                ++Np;
            }
        }
        else
        {
            // If it is a fluid particle calculate based on equation 1 and 2

            // Get an iterator over the neighborhood particles of p
            auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

            // For each neighborhood particle
            while (Np.isNext() == true)
            {
                // ... q
                auto b = Np.get();

                // Get the position xp of the particle
                Point<3, double> xb = vd.getPos(b);

                // if (p == q) skip this particle
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                };

                double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
                Point<3, double> vb = vd.getProp<velocity>(b);
                double Pb = vd.getProp<Pressure>(b);
                double rhob = vd.getProp<rho>(b);

                // Get the distance between p and q
                Point<3, double> dr = xa - xb;
                // take the norm of this vector
                double r2 = norm2(dr);

                // if they interact
                if (r2 < 4.0 * H * H)
                {
                    double r = sqrt(r2);

                    Point<3, double> v_rel = va - vb;

                    Point<3, double> DW;
                    DWab(dr, DW, r, false);

                    double factor = -massb * ((vd.getProp<Pressure>(a) + vd.getProp<Pressure>(b)) / (rhoa * rhob) + Tensile(r, rhoa, rhob, Pa, Pb) + Pi(dr, r2, v_rel, rhoa, rhob, massb, max_visc));

                    vd.getProp<force>(a)[0] += factor * DW.get(0);
                    vd.getProp<force>(a)[1] += factor * DW.get(1);
                    vd.getProp<force>(a)[2] += factor * DW.get(2);

                    vd.getProp<drho>(a) += massb * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
                }

                ++Np;
            }
        }

        ++part;
    }
}

template <typename CellList>
inline void calc_PressureForcesLoop(particles &vd, CellList &NN, double &max_visc)
{
    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();

        // Get the position xp of the particle
        Point<3, double> xa = vd.getPos(a);

        // Take the mass of the particle dependently if it is FLUID or BOUNDARY
        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;

        // Get the density of the of the particle a
        double rhoa = vd.getProp<rho>(a);

        // Get the pressure of the particle a
        double Pa = vd.getProp<Pressure>(a);

        // Get the Velocity of the particle a
        Point<3, double> va = vd.getProp<velocity>(a);

        // Reset the force counter (- gravity on zeta direction)
        vd.template getProp<force>(a)[0] = 0.0;
        vd.template getProp<force>(a)[1] = 0.0;
        vd.template getProp<force>(a)[2] = -gravity;
        vd.template getProp<drho>(a) = 0.0;

        // We threat FLUID particle differently from BOUNDARY PARTICLES ...
        if (vd.getProp<type>(a) != FLUID)
        {
            // If it is a boundary particle calculate the delta rho based on equation 2
            // This require to run across the neighborhoods particles of a
            auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

            // For each neighborhood particle
            while (Np.isNext() == true)
            {
                // ... q
                auto b = Np.get();

                // Get the position xp of the particle
                Point<3, double> xb = vd.getPos(b);

                // if (p == q) skip this particle
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                };

                // get the mass of the particle
                double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;

                // Get the velocity of the particle b
                Point<3, double> vb = vd.getProp<velocity>(b);

                // Get the pressure and density of particle b
                double Pb = vd.getProp<Pressure>(b);
                double rhob = vd.getProp<rho>(b);

                // Get the distance between p and q
                Point<3, double> dr = xa - xb;
                // take the norm of this vector
                double r2 = norm2(dr);

                // If the particles interact ...
                if (r2 < 4.0 * H * H)
                {
                    // ... calculate delta rho
                    double r = sqrt(r2);

                    Point<3, double> dv = va - vb;

                    Point<3, double> DW;
                    DWab(dr, DW, r, false);

                    const double dot = dr.get(0) * dv.get(0) + dr.get(1) * dv.get(1) + dr.get(2) * dv.get(2);
                    const double dot_rr2 = dot / (r2 + Eta2);
                    max_visc = std::max(dot_rr2, max_visc);

                    vd.getProp<drho>(a) += massb * (dv.get(0) * DW.get(0) + dv.get(1) * DW.get(1) + dv.get(2) * DW.get(2));
                }

                ++Np;
            }
        }
        else
        {
            // If it is a fluid particle calculate based on equation 1 and 2

            // Get an iterator over the neighborhood particles of p
            auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

            // For each neighborhood particle
            while (Np.isNext() == true)
            {
                // ... q
                auto b = Np.get();

                // Get the position xp of the particle
                Point<3, double> xb = vd.getPos(b);

                // if (p == q) skip this particle
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                };

                double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
                Point<3, double> vb = vd.getProp<velocity>(b);
                double Pb = vd.getProp<Pressure>(b);
                double rhob = vd.getProp<rho>(b);

                // Get the distance between p and q
                Point<3, double> dr = xa - xb;
                // take the norm of this vector
                double r2 = norm2(dr);

                // if they interact
                if (r2 < 4.0 * H * H)
                {
                    double r = sqrt(r2);

                    Point<3, double> v_rel = va - vb;

                    Point<3, double> DW;
                    DWab(dr, DW, r, false);

                    double factor = -massb * ((vd.getProp<Pressure>(a) + vd.getProp<Pressure>(b)) / (rhoa * rhob));

                    vd.getProp<force>(a)[0] += factor * DW.get(0);
                    vd.getProp<force>(a)[1] += factor * DW.get(1);
                    vd.getProp<force>(a)[2] += factor * DW.get(2);

                    vd.getProp<drho>(a) += massb * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
                }

                ++Np;
            }
        }

        ++part;
    }
}

template <typename CellList>
inline void calc_forces_withoutPressure(particles &vd, CellList &NN, double &max_visc)
{
    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();

        // Get the position xp of the particle
        Point<3, double> xa = vd.getPos(a);

        // Take the mass of the particle dependently if it is FLUID or BOUNDARY
        double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;

        // Get the density of the of the particle a
        double rhoa = vd.getProp<rho>(a);

        // Get the pressure of the particle a
        double Pa = vd.getProp<Pressure>(a);

        // Get the Velocity of the particle a
        Point<3, double> va = vd.getProp<velocity>(a);

        // Reset the force counter (- gravity on zeta direction)
        vd.template getProp<force>(a)[0] = 0.0;
        vd.template getProp<force>(a)[1] = 0.0;
        vd.template getProp<force>(a)[2] = -gravity;
        vd.template getProp<drho>(a) = 0.0;

        // We threat FLUID particle differently from BOUNDARY PARTICLES ...
        if (vd.getProp<type>(a) != FLUID)
        {
            // If it is a boundary particle calculate the delta rho based on equation 2
            // This require to run across the neighborhoods particles of a

            auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

            // For each neighborhood particle
            while (Np.isNext() == true)
            {
                // ... q
                auto b = Np.get();

                // Get the position xp of the particle
                Point<3, double> xb = vd.getPos(b);

                // if (p == q) skip this particle
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                };

                // get the mass of the particle
                double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;

                // Get the velocity of the particle b
                Point<3, double> vb = vd.getProp<velocity>(b);

                // Get the pressure and density of particle b
                double Pb = vd.getProp<Pressure>(b);
                double rhob = vd.getProp<rho>(b);

                // Get the distance between p and q
                Point<3, double> dr = xa - xb;
                // take the norm of this vector
                double r2 = norm2(dr);

                // If the particles interact ...
                if (r2 < 4.0 * H * H)
                {
                    // ... calculate delta rho
                    double r = sqrt(r2);

                    Point<3, double> dv = va - vb;

                    Point<3, double> DW;
                    DWab(dr, DW, r, false);

                    const double dot = dr.get(0) * dv.get(0) + dr.get(1) * dv.get(1) + dr.get(2) * dv.get(2);
                    const double dot_rr2 = dot / (r2 + Eta2);
                    max_visc = std::max(dot_rr2, max_visc);

                    vd.getProp<drho>(a) += massb * (dv.get(0) * DW.get(0) + dv.get(1) * DW.get(1) + dv.get(2) * DW.get(2));
                }

                ++Np;
            }
        }
        else
        {
            // If it is a fluid particle calculate based on equation 1 and 2

            // Get an iterator over the neighborhood particles of p
            auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

            // For each neighborhood particle
            while (Np.isNext() == true)
            {
                // ... q
                auto b = Np.get();

                // Get the position xp of the particle
                Point<3, double> xb = vd.getPos(b);

                // if (p == q) skip this particle
                if (a.getKey() == b)
                {
                    ++Np;
                    continue;
                };

                double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
                Point<3, double> vb = vd.getProp<velocity>(b);
                double Pb = vd.getProp<Pressure>(b);
                double rhob = vd.getProp<rho>(b);

                // Get the distance between p and q
                Point<3, double> dr = xa - xb;
                // take the norm of this vector
                double r2 = norm2(dr);

                // if they interact
                if (r2 < 4.0 * H * H)
                {
                    double r = sqrt(r2);

                    Point<3, double> v_rel = va - vb;

                    Point<3, double> DW;
                    DWab(dr, DW, r, false);

                    // double factor = -massb * ((vd.getProp<Pressure>(a) + vd.getProp<Pressure>(b)) / (rhoa * rhob) + Tensile(r, rhoa, rhob, Pa, Pb) + Pi(dr, r2, v_rel, rhoa, rhob, massb, max_visc));
                    double factor = -massb * (+Tensile(r, rhoa, rhob, Pa, Pb) + Pi(dr, r2, v_rel, rhoa, rhob, massb, max_visc));

                    vd.getProp<force>(a)[0] += factor * DW.get(0);
                    vd.getProp<force>(a)[1] += factor * DW.get(1);
                    vd.getProp<force>(a)[2] += factor * DW.get(2);

                    vd.getProp<drho>(a) += massb * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
                }

                ++Np;
            }
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

        Point<3, double> acc(vd.getProp<force>(a));
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
    double dt = double(CFLnumber) * std::min(dt_f, dt_cv);
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
            vd.template getProp<velocity>(a)[0] = 0.0;
            vd.template getProp<velocity>(a)[1] = 0.0;
            vd.template getProp<velocity>(a)[2] = 0.0;
            double rhonew = vd.template getProp<rho_prev>(a) + dt2 * vd.template getProp<drho>(a);
            vd.template getProp<rho>(a) = (rhonew < rho_zero) ? rho_zero : rhonew;

            vd.template getProp<rho_prev>(a) = rhop;

            ++part;
            continue;
        }

        //-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
        double dx = vd.template getProp<velocity>(a)[0] * dt + vd.template getProp<force>(a)[0] * dt205;
        double dy = vd.template getProp<velocity>(a)[1] * dt + vd.template getProp<force>(a)[1] * dt205;
        double dz = vd.template getProp<velocity>(a)[2] * dt + vd.template getProp<force>(a)[2] * dt205;

        vd.getPos(a)[0] += dx;
        vd.getPos(a)[1] += dy;
        vd.getPos(a)[2] += dz;

        double velX = vd.template getProp<velocity>(a)[0];
        double velY = vd.template getProp<velocity>(a)[1];
        double velZ = vd.template getProp<velocity>(a)[2];
        double rhop = vd.template getProp<rho>(a);

        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + vd.template getProp<force>(a)[0] * dt2;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + vd.template getProp<force>(a)[1] * dt2;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + vd.template getProp<force>(a)[2] * dt2;
        vd.template getProp<rho>(a) = vd.template getProp<rho_prev>(a) + dt2 * vd.template getProp<drho>(a);

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
            // Update rho
            double rhop = vd.template getProp<rho>(a);

            // Update only the density
            vd.template getProp<velocity>(a)[0] = 0.0;
            vd.template getProp<velocity>(a)[1] = 0.0;
            vd.template getProp<velocity>(a)[2] = 0.0;
            double rhonew = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);
            vd.template getProp<rho>(a) = (rhonew < rho_zero) ? rho_zero : rhonew;
            vd.template getProp<rho_prev>(a) = rhop;
            ++part;
            continue;
        }
        //-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
        double dx = vd.template getProp<velocity>(a)[0] * dt + vd.template getProp<force>(a)[0] * dt205;
        double dy = vd.template getProp<velocity>(a)[1] * dt + vd.template getProp<force>(a)[1] * dt205;
        double dz = vd.template getProp<velocity>(a)[2] * dt + vd.template getProp<force>(a)[2] * dt205;
        vd.getPos(a)[0] += dx;
        vd.getPos(a)[1] += dy;
        vd.getPos(a)[2] += dz;
        double velX = vd.template getProp<velocity>(a)[0];
        double velY = vd.template getProp<velocity>(a)[1];
        double velZ = vd.template getProp<velocity>(a)[2];
        double rhop = vd.template getProp<rho>(a);
        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity>(a)[0] + vd.template getProp<force>(a)[0] * dt;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity>(a)[1] + vd.template getProp<force>(a)[1] * dt;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity>(a)[2] + vd.template getProp<force>(a)[2] * dt;
        vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);
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

void verlet_int_no_densityUpdate(particles &vd, double dt)
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
            vd.template getProp<velocity>(a)[0] = 0.0;
            vd.template getProp<velocity>(a)[1] = 0.0;
            vd.template getProp<velocity>(a)[2] = 0.0;
            /*

            double rhonew = vd.template getProp<rho_prev>(a) + dt2 * vd.template getProp<drho>(a);
            vd.template getProp<rho>(a) = (rhonew < rho_zero) ? rho_zero : rhonew;

            */
            vd.template getProp<rho_prev>(a) = rhop;

            ++part;
            continue;
        }

        //-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
        double dx = vd.template getProp<velocity>(a)[0] * dt + vd.template getProp<force>(a)[0] * dt205;
        double dy = vd.template getProp<velocity>(a)[1] * dt + vd.template getProp<force>(a)[1] * dt205;
        double dz = vd.template getProp<velocity>(a)[2] * dt + vd.template getProp<force>(a)[2] * dt205;

        vd.getPos(a)[0] += dx;
        vd.getPos(a)[1] += dy;
        vd.getPos(a)[2] += dz;

        double velX = vd.template getProp<velocity>(a)[0];
        double velY = vd.template getProp<velocity>(a)[1];
        double velZ = vd.template getProp<velocity>(a)[2];
        double rhop = vd.template getProp<rho>(a);

        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + vd.template getProp<force>(a)[0] * dt2;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + vd.template getProp<force>(a)[1] * dt2;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + vd.template getProp<force>(a)[2] * dt2;
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

void euler_intno_densityUpdate(particles &vd, double dt)
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
            vd.template getProp<velocity>(a)[0] = 0.0;
            vd.template getProp<velocity>(a)[1] = 0.0;
            vd.template getProp<velocity>(a)[2] = 0.0;
            /*
            double rhonew = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);
            vd.template getProp<rho>(a) = (rhonew < rho_zero) ? rho_zero : rhonew;
            */
            vd.template getProp<rho_prev>(a) = rhop;
            ++part;
            continue;
        }
        //-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
        double dx = vd.template getProp<velocity>(a)[0] * dt + vd.template getProp<force>(a)[0] * dt205;
        double dy = vd.template getProp<velocity>(a)[1] * dt + vd.template getProp<force>(a)[1] * dt205;
        double dz = vd.template getProp<velocity>(a)[2] * dt + vd.template getProp<force>(a)[2] * dt205;
        vd.getPos(a)[0] += dx;
        vd.getPos(a)[1] += dy;
        vd.getPos(a)[2] += dz;
        double velX = vd.template getProp<velocity>(a)[0];
        double velY = vd.template getProp<velocity>(a)[1];
        double velZ = vd.template getProp<velocity>(a)[2];
        double rhop = vd.template getProp<rho>(a);
        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity>(a)[0] + vd.template getProp<force>(a)[0] * dt;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity>(a)[1] + vd.template getProp<force>(a)[1] * dt;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity>(a)[2] + vd.template getProp<force>(a)[2] * dt;
        // vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);
        //  Check if the particle go out of range in space and in density
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
inline void EqState_incompressible(particles &vd, CellList &NN, double &dt, double &max_visc)
{
    double error = 0;
    double error_min = 0.01;
    int min_itteration = 3;
    int iter = 0;

    double dt205 = dt * dt * 0.5;
    double dt2 = dt * 2.0;
    double rho_predicted_variation;

    auto it = vd.getDomainIterator();

    // whiler max_error > something do
    while ((error > error_min) || (iter < min_itteration))
    {
        while (it.isNext())
        {
            /* Velocity and Position block */
            auto a = it.get();

            Point<3, double> xa = vd.getPos(a);
            // Take the mass of the particle dependently if it is FLUID or BOUNDARY
            double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;

            // Get the density of the of the particle a
            double rhoa = vd.getProp<rho>(a);

            // Get the pressure of the particle a
            double Pa = vd.getProp<Pressure>(a);

            // Get the Velocity of the particle a
            Point<3, double> va = vd.getProp<velocity>(a);

            double fx_without_pressure = vd.template getProp<force>(a)[0];
            double fy_without_pressure = vd.template getProp<force>(a)[1];
            double fz_without_pressure = vd.template getProp<force>(a)[2];

            double fx_pressure_force = 0;
            double fy_pressure_force = 0;
            double fz_pressure_force = 0;

            double pressure = 0;
            vd.template getProp<drho>(a) = 0.0;

            // 	get and update velocity

            vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity>(a)[0] + (fx_without_pressure + fx_pressure_force) * dt;
            vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity>(a)[1] + (fy_without_pressure + fy_pressure_force) * dt;
            vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity>(a)[2] + (fz_without_pressure + fz_pressure_force) * dt;

            //	Calculate displacement and update position

            double dx = vd.template getProp<velocity>(a)[0] * dt + (fx_without_pressure + fx_pressure_force) * dt205;
            double dy = vd.template getProp<velocity>(a)[1] * dt + (fy_without_pressure + fy_pressure_force) * dt205;
            double dz = vd.template getProp<velocity>(a)[2] * dt + (fz_without_pressure + fz_pressure_force) * dt205;

            vd.getPos(a)[0] += dx;
            vd.getPos(a)[1] += dy;
            vd.getPos(a)[2] += dz;

            // Check if the particle go out of range in space and in density
            /*
            need to check if this is necesarry
            */
            if (vd.getPos(a)[0] < 0.000263878 || vd.getPos(a)[1] < 0.000263878 || vd.getPos(a)[2] < 0.000263878 ||
                vd.getPos(a)[0] > 0.000263878 + 1.59947 || vd.getPos(a)[1] > 0.000263878 + 0.672972 || vd.getPos(a)[2] > 0.000263878 + 0.903944 ||
                vd.template getProp<rho>(a) < RhoMin || vd.template getProp<rho>(a) > RhoMax)
            {
                to_remove.add(a.getKey());
            }
            /* Density update */

            if (vd.getProp<type>(a) != FLUID)
            {
                // If it is a boundary particle calculate the delta rho based on equation 2
                // This require to run across the neighborhoods particles of a
                auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

                // For each neighborhood particle
                while (Np.isNext() == true)
                {
                    // ... q
                    auto b = Np.get();

                    // Get the position xp of the particle
                    Point<3, double> xb = vd.getPos(b);

                    // if (p == q) skip this particle
                    if (a.getKey() == b)
                    {
                        ++Np;
                        continue;
                    };

                    // get the mass of the particle
                    double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;

                    // Get the velocity of the particle b
                    Point<3, double> vb = vd.getProp<velocity>(b);

                    // Get the pressure and density of particle b
                    double Pb = vd.getProp<Pressure>(b);
                    double rhob = vd.getProp<rho>(b);

                    // Get the distance between p and q
                    Point<3, double> dr = xa - xb;

                    // take the norm of this vector
                    double r2 = norm2(dr);

                    // If the particles interact ...
                    if (r2 < 4.0 * H * H)
                    {
                        // ... calculate delta rho
                        double r = sqrt(r2);
                        Wab(r);
                        /*
                        Point<3, double> dv = va - vb;

                        Point<3, double> DW;

                        DWab(dr, DW, r, false);

                        const double dot = dr.get(0) * dv.get(0) + dr.get(1) * dv.get(1) + dr.get(2) * dv.get(2);
                        const double dot_rr2 = dot / (r2 + Eta2);
                        max_visc = std::max(dot_rr2, max_visc);

                        vd.getProp<drho>(a) += massb * (dv.get(0) * DW.get(0) + dv.get(1) * DW.get(1) + dv.get(2) * DW.get(2));
                        */
                        vd.getProp<rho>(a) += massb * (Wab(r));
                    }
                    ++Np;
                }
            }
            else
            {
                // Get an iterator over the neighborhood particles of p
                auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

                // For each neighborhood particle
                while (Np.isNext() == true)
                {
                    // ... q
                    auto b = Np.get();

                    // Get the position xp of the particle
                    Point<3, double> xb = vd.getPos(b);

                    // if (p == q) skip this particle
                    if (a.getKey() == b)
                    {
                        ++Np;
                        continue;
                    };

                    double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;
                    Point<3, double> vb = vd.getProp<velocity>(b);
                    double Pb = vd.getProp<Pressure>(b);
                    double rhob = vd.getProp<rho>(b);

                    // Get the distance between p and q
                    Point<3, double> dr = xa - xb;
                    // take the norm of this vector
                    double r2 = norm2(dr);

                    // if they interact
                    if (r2 < 4.0 * H * H)
                    {
                        double r = sqrt(r2);

                        Point<3, double> v_rel = va - vb;

                        Point<3, double> DW;
                        DWab(dr, DW, r, false);

                        vd.getProp<drho>(a) += massb * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
                    }
                    ++Np;
                }
            }
            double rho_predicted = vd.getProp<rho>(a) + vd.getProp<drho>(a) * dt2;
            rho_predicted_variation = std::abs((rho_predicted - rho_zero) / rho_zero);
        };
        error = std::max(error, rho_predicted_variation);
    };
};
template <typename CellList>
inline void calc_Density(particles &vd, CellList &NN)
{
    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto a = it.get();
        //  initialisiere mit 0
        vd.getProp<rho>(a) = 0;
        Point<3, double> xa = vd.getPos(a);
        // Take the mass of the particle dependently if it is FLUID or BOUNDARY
        double massa = MassFluid;
        /* Density update */
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
        // For each neighborhood particle
        while (Np.isNext() == true)
        {
            auto b = Np.get();
            Point<3, double> xb = vd.getPos(b);
            if (a.getKey() == b)
            {
                ++Np;
                continue;
            };
            double r = sqrt(norm2(xa - xb));
            LOGFunction("r : ");
            LOGFunction(r);

            LOGFunction("calling Wab : ");

            LOGFunction(Wab(r));
            double ker = Wab(r)     ; 

            vd.getProp<rho>(a) += ker;
            ++Np;
        }

        vd.getProp<rho>(a) = vd.getProp<rho>(a) * massa;
        LOGFunction("rho : ");
        LOGFunction(vd.getProp<rho>(a));
        if (vd.template getProp<rho>(a) < RhoMin || vd.template getProp<rho>(a) > RhoMax)
        {
            to_remove.add(a.getKey());
        }
        ++it;
    };
};