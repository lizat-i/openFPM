#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"

typedef vector_dist<3, double, aggregate<size_t, double, double, double, double, double[3], double[3], double[3]>> particles;

const size_t type = 0;
struct operationParams
{
    const int BOUNDARY;
    const int FLUID;
    const double dp;
    double h_swl;
    const double coeff_sound;
    const double gamma_;
    const double H;
    const double Eta2;
    const double visco;
    double cbar;
    const double MassFluid;
    const double MassBound;
    const double gravity;
    const double t_end;
    const double rho_zero;
    double B;
    double W_dap;
    const double CFLnumber;
    const double DtMin;
    const double RhoMin;
    const double RhoMax;
    double max_fluid_height;
    const double a2;
    const double c1;
    const double d1;
    const double c2;
    const double a2_4;
    size_t cnt;
    double max_visc;
};
enum properties
{
  rho = 1,
  rho_prev,
  Pressure,
  drho,
  force,
  velocity,
  velocity_prev
};
openfpm::vector<size_t> to_remove;

inline void EqState(particles & vd, operationParams &params)
{    
    auto it = vd.getDomainIterator();

    while (it.isNext())
    {
        auto a = it.get();
        double rho_a = vd.template getProp<rho>(a);
        double rho_frac = rho_a / params.rho_zero;
        vd.template getProp<Pressure>(a) = params.B*( rho_frac*rho_frac*rho_frac*rho_frac*rho_frac*rho_frac*rho_frac - 1.0);
        ++it;
    }
}
inline double Wab(double r, operationParams &params)
{
    r /= params.H;
    if (r < 1.0)
        return (1.0 - 3.0/2.0*r*r + 3.0/4.0*r*r*r)*params.a2;
    else if (r < 2.0)
        return (1.0/4.0*(2.0 - r*r)*(2.0 - r*r)*(2.0 - r*r))*params.a2;
    else
        return 0.0;
}

inline void DWab(Point<3,double> & dx, Point<3,double> & DW, double r, bool print, operationParams &params )
{   

    const double qq=r/params.H;
    double qq2 = qq * qq;
    double fac1 = (params.c1*qq + params.d1*qq2)/r;
    double b1 = (qq < 1.0)?1.0f:0.0f;
    double wqq = (2.0 - qq);
    double fac2 = params.c2 * wqq * wqq / r;
    double b2 = (qq >= 1.0 && qq < 2.0)?1.0f:0.0f;
    double factor = (b1*fac1 + b2*fac2);
    DW.get(0) = factor * dx.get(0);
    DW.get(1) = factor * dx.get(1);
    DW.get(2) = factor * dx.get(2);
}
// Tensile correction
inline double Tensile(double r, double rhoa, double rhob, double prs1, double prs2, operationParams &params)
{
    const double qq=r/params.H; 
    double wab;

    if(r>params.H)
    {
        double wqq1=2.0f-qq;
        double wqq2=wqq1*wqq1;
        wab=params.a2_4*(wqq2*wqq1);
    }
    else
    {
        double wqq2=qq*qq;
        double wqq3=wqq2*qq;
        wab=params.a2*(1.0f-1.5f*wqq2+0.75f*wqq3);
    }
    //-Tensile correction.
    double fab=wab*params.W_dap;
    fab*=fab; fab*=fab; //fab=fab^4
    const double tensilp1=(prs1/(rhoa*rhoa))*(prs1>0? 0.01: -0.2);
    const double tensilp2=(prs2/(rhob*rhob))*(prs2>0? 0.01: -0.2);
    return (fab*(tensilp1+tensilp2));
}
inline double Pi(const Point<3,double> & dr, double rr2, Point<3,double> & dv, double rhoa, double rhob, double massb, double & visc, operationParams &params)
{
    const double dot = dr.get(0)*dv.get(0) + dr.get(1)*dv.get(1) + dr.get(2)*dv.get(2);
    const double dot_rr2 = dot/(rr2+params.Eta2);
    visc=std::max(dot_rr2,visc);
    if(dot < 0)
    {
        const float amubar=params.H*dot_rr2;
        const float robar=(rhoa+rhob)*0.5f;
        const float pi_visc=(-params.visco*params.cbar*amubar/robar);
        return pi_visc;
    }
    else
        return 0.0;
}

template<typename CellList> inline void calc_forces(particles & vd, CellList & NN, double & max_visc, operationParams &params)
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
        Point<3,double> xa = vd.getPos(a);
        // Take the mass of the particle dependently if it is FLUID or BOUNDARY
        double massa = (vd.getProp<type>(a) == params.FLUID)?params.MassFluid:params.MassBound;
        // Get the density of the of the particle a
        double rhoa = vd.getProp<rho>(a);
        // Get the pressure of the particle a
        double Pa = vd.getProp<Pressure>(a);
        // Get the Velocity of the particle a
        Point<3,double> va = vd.getProp<velocity>(a);
        // Reset the force counter (- gravity on zeta direction)
        vd.template getProp<force>(a)[0] = 0.0;
        vd.template getProp<force>(a)[1] = 0.0;
        vd.template getProp<force>(a)[2] = -params.gravity;
        vd.template getProp<drho>(a) = 0.0;
        // We threat FLUID particle differently from BOUNDARY PARTICLES ...
        if (vd.getProp<type>(a) != params.FLUID)
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
                Point<3,double> xb = vd.getPos(b);
                // if (p == q) skip this particle
                if (a.getKey() == b)    {++Np; continue;};
                // get the mass of the particle
                double massb = (vd.getProp<type>(b) == params.FLUID)?params.MassFluid:params.MassBound;
                // Get the velocity of the particle b
                Point<3,double> vb = vd.getProp<velocity>(b);
                // Get the pressure and density of particle b
                double Pb = vd.getProp<Pressure>(b);
                double rhob = vd.getProp<rho>(b);
                // Get the distance between p and q
                Point<3,double> dr = xa - xb;
                // take the norm of this vector
                double r2 = norm2(dr);
                // If the particles interact ...
                if (r2 < 4.0*params.H*params.a2)
                {
                    // ... calculate delta rho
                    double r = sqrt(r2);
                    Point<3,double> dv = va - vb;
                    Point<3,double> DW;
                    DWab(dr,DW,r,false,params);
                    const double dot = dr.get(0)*dv.get(0) + dr.get(1)*dv.get(1) + dr.get(2)*dv.get(2);
                    const double dot_rr2 = dot/(r2+params.Eta2);
                    max_visc=std::max(dot_rr2,params.max_visc);
                    vd.getProp<drho>(a) += massb*(dv.get(0)*DW.get(0)+dv.get(1)*DW.get(1)+dv.get(2)*DW.get(2));
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
                Point<3,double> xb = vd.getPos(b);
                // if (p == q) skip this particle
                if (a.getKey() == b)    {++Np; continue;};
                double massb = (vd.getProp<type>(b) == params.FLUID)?params.MassFluid:params.MassBound;
                Point<3,double> vb = vd.getProp<velocity>(b);
                double Pb = vd.getProp<Pressure>(b);
                double rhob = vd.getProp<rho>(b);
                // Get the distance between p and q
                Point<3,double> dr = xa - xb;
                // take the norm of this vector
                double r2 = norm2(dr);
                // if they interact
                if (r2 < 4.0*params.H*params.H)
                {
                    double r = sqrt(r2);
                    Point<3,double> v_rel = va - vb;
                    Point<3,double> DW;
                    DWab(dr,DW,r,false,params);
                    double factor = - massb*((vd.getProp<Pressure>(a) + vd.getProp<Pressure>(b)) / (rhoa * rhob) + Tensile(r,rhoa,rhob,Pa,Pb,params) + Pi(dr,r2,v_rel,rhoa,rhob,massb,max_visc,params));
                    vd.getProp<force>(a)[0] += factor * DW.get(0);
                    vd.getProp<force>(a)[1] += factor * DW.get(1);
                    vd.getProp<force>(a)[2] += factor * DW.get(2);
                    vd.getProp<drho>(a) += massb*(v_rel.get(0)*DW.get(0)+v_rel.get(1)*DW.get(1)+v_rel.get(2)*DW.get(2));
                }
                ++Np;
            }
        }
        ++part;
    }
}

void max_acceleration_and_velocity(particles & vd, double & max_acc, double & max_vel)
{
    // Calculate the maximum acceleration
    auto part = vd.getDomainIterator();
    while (part.isNext())
    {
        auto a = part.get();
        Point<3,double> acc(vd.getProp<force>(a));
        double acc2 = norm2(acc);
        Point<3,double> vel(vd.getProp<velocity>(a));
        double vel2 = norm2(vel);
        if (vel2 >= max_vel)
            max_vel = vel2;
        if (acc2 >= max_acc)
            max_acc = acc2;
        ++part;
    }
    max_acc = sqrt(max_acc);
    max_vel = sqrt(max_vel);
    Vcluster<> & v_cl = create_vcluster();
    v_cl.max(max_acc);
    v_cl.max(max_vel);
    v_cl.execute();
}

double calc_deltaT(particles & vd, double ViscDtMax,  operationParams &params)
{
    double Maxacc = 0.0;
    double Maxvel = 0.0;
    max_acceleration_and_velocity(vd,Maxacc,Maxvel);
    //-dt1 depends on force per unit mass.
    const double dt_f = (Maxacc)?sqrt(params.H/Maxacc):std::numeric_limits<int>::max();
    //-dt2 combines the Courant and the viscous time-step controls.
    const double dt_cv = params.H/(std::max(params.cbar,Maxvel*10.) + params.H*ViscDtMax);
    //-dt new value of time step.
    double dt=double(params.CFLnumber)*std::min(dt_f,dt_cv);
    if(dt<double(params.DtMin))
        dt=double(params.DtMin);
    return dt;
}
void verlet_int(particles & vd, double dt, operationParams &params)
{
    // list of the particle to remove
    to_remove.clear();
    // particle iterator
    auto part = vd.getDomainIterator();
    double dt205 = dt*dt*0.5;
    double dt2 = dt*2.0;
    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();
        // if the particle is boundary
        if (vd.template getProp<type>(a) == params.BOUNDARY)
        {
            // Update rho
            double rhop = vd.template getProp<rho>(a);
            // Update only the density
            vd.template getProp<velocity>(a)[0] = 0.0;
            vd.template getProp<velocity>(a)[1] = 0.0;
            vd.template getProp<velocity>(a)[2] = 0.0;
            double rhonew = vd.template getProp<rho_prev>(a) + dt2*vd.template getProp<drho>(a);
            vd.template getProp<rho>(a) = (rhonew < params.rho_zero)?params.rho_zero:rhonew;
            vd.template getProp<rho_prev>(a) = rhop;
            ++part;
            continue;
        }
        //-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
        double dx = vd.template getProp<velocity>(a)[0]*dt + vd.template getProp<force>(a)[0]*dt205;
        double dy = vd.template getProp<velocity>(a)[1]*dt + vd.template getProp<force>(a)[1]*dt205;
        double dz = vd.template getProp<velocity>(a)[2]*dt + vd.template getProp<force>(a)[2]*dt205;
        vd.getPos(a)[0] += dx;
        vd.getPos(a)[1] += dy;
        vd.getPos(a)[2] += dz;
        double velX = vd.template getProp<velocity>(a)[0];
        double velY = vd.template getProp<velocity>(a)[1];
        double velZ = vd.template getProp<velocity>(a)[2];
        double rhop = vd.template getProp<rho>(a);
        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + vd.template getProp<force>(a)[0]*dt2;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + vd.template getProp<force>(a)[1]*dt2;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + vd.template getProp<force>(a)[2]*dt2;
        vd.template getProp<rho>(a) = vd.template getProp<rho_prev>(a) + dt2*vd.template getProp<drho>(a);
        // Check if the particle go out of range in space and in density
        if (vd.getPos(a)[0] <  0.000263878 || vd.getPos(a)[1] < 0.000263878 || vd.getPos(a)[2] < 0.000263878 ||
            vd.getPos(a)[0] >  0.000263878+1.59947 || vd.getPos(a)[1] > 0.000263878+0.672972 || vd.getPos(a)[2] > 0.000263878+0.903944 ||
            vd.template getProp<rho>(a) < params.RhoMin || vd.template getProp<rho>(a) > params.RhoMax)
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
    vd.remove(to_remove,0);
    // increment the iteration counter
    params.cnt++;
}
void euler_int(particles & vd, double dt,operationParams &params )
{
    // list of the particle to remove
    to_remove.clear();
    // particle iterator
    auto part = vd.getDomainIterator();
    double dt205 = dt*dt*0.5;
    double dt2 = dt*2.0;
    // For each particle ...
    while (part.isNext())
    {
        // ... a
        auto a = part.get();
        // if the particle is boundary
        if (vd.template getProp<type>(a) == params.BOUNDARY)
        {
            // Update rho
            double rhop = vd.template getProp<rho>(a);
            // Update only the density
            vd.template getProp<velocity>(a)[0] = 0.0;
            vd.template getProp<velocity>(a)[1] = 0.0;
            vd.template getProp<velocity>(a)[2] = 0.0;
            double rhonew = vd.template getProp<rho>(a) + dt*vd.template getProp<drho>(a);
            vd.template getProp<rho>(a) = (rhonew < params.rho_zero)?params.rho_zero:rhonew;
            vd.template getProp<rho_prev>(a) = rhop;
            ++part;
            continue;
        }
        //-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
        double dx = vd.template getProp<velocity>(a)[0]*dt + vd.template getProp<force>(a)[0]*dt205;
        double dy = vd.template getProp<velocity>(a)[1]*dt + vd.template getProp<force>(a)[1]*dt205;
        double dz = vd.template getProp<velocity>(a)[2]*dt + vd.template getProp<force>(a)[2]*dt205;
        vd.getPos(a)[0] += dx;
        vd.getPos(a)[1] += dy;
        vd.getPos(a)[2] += dz;
        double velX = vd.template getProp<velocity>(a)[0];
        double velY = vd.template getProp<velocity>(a)[1];
        double velZ = vd.template getProp<velocity>(a)[2];
        double rhop = vd.template getProp<rho>(a);
        vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity>(a)[0] + vd.template getProp<force>(a)[0]*dt;
        vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity>(a)[1] + vd.template getProp<force>(a)[1]*dt;
        vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity>(a)[2] + vd.template getProp<force>(a)[2]*dt;
        vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt*vd.template getProp<drho>(a);
        // Check if the particle go out of range in space and in density
        if (vd.getPos(a)[0] <  0.000263878 || vd.getPos(a)[1] < 0.000263878 || vd.getPos(a)[2] < 0.000263878 ||
            vd.getPos(a)[0] >  0.000263878+1.59947 || vd.getPos(a)[1] > 0.000263878+0.672972 || vd.getPos(a)[2] > 0.000263878+0.903944 ||
            vd.template getProp<rho>(a) < params.RhoMin || vd.template getProp<rho>(a) > params.RhoMax)
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
    vd.remove(to_remove,0);
    // increment the iteration counter
    params.cnt++;
}

template<typename Vector, typename CellList>
inline void sensor_pressure(Vector & vd,
                            CellList & NN,
                            openfpm::vector<openfpm::vector<double>> & press_t,
                            openfpm::vector<Point<3,double>> & probes,
                            operationParams &params)
{
    Vcluster<> & v_cl = create_vcluster();
    press_t.add();
    for (size_t i = 0 ; i < probes.size() ; i++)
    {
        float press_tmp = 0.0f;
        float tot_ker = 0.0;
        // if the probe is inside the processor domain
        if (vd.getDecomposition().isLocal(probes.get(i)) == true)
        {
            // Get the position of the probe i
            Point<3,double> xp = probes.get(i);
            // get the iterator over the neighbohood particles of the probes position
            auto itg = NN.template getNNIterator<NO_CHECK>(NN.getCell(probes.get(i)));
            while (itg.isNext())
            {
                auto q = itg.get();
                // Only the fluid particles are importants
                if (vd.template getProp<type>(q) != params.FLUID)
                {
                    ++itg;
                    continue;
                }
                // Get the position of the neighborhood particle q
                Point<3,double> xq = vd.getPos(q);
                // Calculate the contribution of the particle to the pressure
                // of the probe
                double r = sqrt(norm2(xp - xq));
                double ker = Wab(r,params) * (params.MassFluid / params.rho_zero);
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
struct ModelCustom
{
    template <typename Decomposition, typename vector>
    inline void addComputation(Decomposition &dec,
                               vector &vd,
                               size_t v,
                               size_t p)
    {
        if (vd.template getProp<type>(p) == 1)
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

void createBoxAndParseBox(particles &vd, operationParams &parameters, Box<3,double> &domain, Ghost<3,double> &g)
{   size_t sz[3] = {207, 90, 66};
    Box<3, double> fluid_box({parameters.dp / 2.0, parameters.dp / 2.0, parameters.dp / 2.0}, {0.4 + parameters.dp / 2.0, 0.67 - parameters.dp / 2.0, 0.3 + parameters.dp / 2.0});
    auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);
    // here we fill some of the constants needed by the simulation
    parameters.max_fluid_height = fluid_it.getBoxMargins().getHigh(2);
    parameters.h_swl = fluid_it.getBoxMargins().getHigh(2) - fluid_it.getBoxMargins().getLow(2);
    parameters.B = (parameters.coeff_sound) * (parameters.coeff_sound) * parameters.gravity * parameters.h_swl * parameters.rho_zero / parameters.gamma_;
    parameters.cbar = parameters.coeff_sound * sqrt(parameters.gravity * parameters.h_swl);
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
        vd.template getLastProp<type>() = parameters.FLUID;
        // We also initialize the density of the particle and the hydro-static pressure given by
        //
        // rho_zero*g*h = P
        //
        // rho_p = (P/B + 1)^(1/Gamma) * rho_zero
        //
        vd.template getLastProp<Pressure>() = parameters.rho_zero * parameters.gravity * (parameters.max_fluid_height - fluid_it.get().get(2));
        vd.template getLastProp<rho>() = pow(vd.template getLastProp<Pressure>() / parameters.B + 1, 1.0 / parameters.gamma_) * parameters.rho_zero;
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
    Box<3, double> recipient1({0.0, 0.0, 0.0}, {1.6 + parameters.dp / 2.0, 0.67 + parameters.dp / 2.0, 0.4 + parameters.dp / 2.0});
    Box<3, double> recipient2({parameters.dp, parameters.dp, parameters.dp}, {1.6 - parameters.dp / 2.0, 0.67 - parameters.dp / 2.0, 0.4 + parameters.dp / 2.0});
    Box<3, double> obstacle1({0.9, 0.24 - parameters.dp / 2.0, 0.0}, {1.02 + parameters.dp / 2.0, 0.36, 0.45 + parameters.dp / 2.0});
    Box<3, double> obstacle2({0.9 + parameters.dp, 0.24 + parameters.dp / 2.0, 0.0}, {1.02 - parameters.dp / 2.0, 0.36 - parameters.dp, 0.45 - parameters.dp / 2.0});
    Box<3, double> obstacle3({0.9 + parameters.dp, 0.24, 0.0}, {1.02, 0.36, 0.45});
    openfpm::vector<Box<3, double>> holes;
    holes.add(recipient2);
    holes.add(obstacle1);
    auto bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient1);
    while (bound_box.isNext())
    {
        vd.add();
        vd.getLastPos()[0] = bound_box.get().get(0);
        vd.getLastPos()[1] = bound_box.get().get(1);
        vd.getLastPos()[2] = bound_box.get().get(2);
        vd.template getLastProp<type>() = parameters.BOUNDARY;
        vd.template getLastProp<rho>() = parameters.rho_zero;
        vd.template getLastProp<rho_prev>() = parameters.rho_zero;
        vd.template getLastProp<velocity>()[0] = 0.0;
        vd.template getLastProp<velocity>()[1] = 0.0;
        vd.template getLastProp<velocity>()[2] = 0.0;
        vd.template getLastProp<velocity_prev>()[0] = 0.0;
        vd.template getLastProp<velocity_prev>()[1] = 0.0;
        vd.template getLastProp<velocity_prev>()[2] = 0.0;
        ++bound_box;
    }
    auto obstacle_box = DrawParticles::DrawSkin(vd, sz, domain, obstacle2, obstacle1);
    while (obstacle_box.isNext())
    {
        vd.add();
        vd.getLastPos()[0] = obstacle_box.get().get(0);
        vd.getLastPos()[1] = obstacle_box.get().get(1);
        vd.getLastPos()[2] = obstacle_box.get().get(2);
        vd.template getLastProp<type>() = parameters.BOUNDARY;
        vd.template getLastProp<rho>() = parameters.rho_zero;
        vd.template getLastProp<rho_prev>() = parameters.rho_zero;
        vd.template getLastProp<velocity>()[0] = 0.0;
        vd.template getLastProp<velocity>()[1] = 0.0;
        vd.template getLastProp<velocity>()[2] = 0.0;
        vd.template getLastProp<velocity_prev>()[0] = 0.0;
        vd.template getLastProp<velocity_prev>()[1] = 0.0;
        vd.template getLastProp<velocity_prev>()[2] = 0.0;
        ++obstacle_box;
    }
}