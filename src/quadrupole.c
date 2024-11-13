#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

static void rebx_calculate_quad_force(struct reb_simulation* const sim, struct reb_particle* const particles, const int N, const double mu_eff, const double R, const int source_index){
    const struct reb_particle source = particles[source_index];
//    can change to hardcode?
    for (int i=0; i<N; i++){
        if(i != 0){
            continue;
        }
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        
//      Hard coding gamma and A
        const double gamma = -4.0;
        const double G = sim->G;
        const double A = 3 * G * p.m / (4 * source.m);
        const double prefac = mu_eff*A*pow(R,2)*pow(r2, (gamma-1.)/2.);

        // reversing the direction of force: apply to "source" (Earth) FROM particle "i" (Sun)
        particles[source_index].ax += prefac * dx;
        particles[source_index].ay += prefac * dy;
        particles[source_index].az += prefac * dz;

        // apply opposite force to particle "i" (Sun)
        particles[i].ax -= source.m/p.m * prefac * dx;
        particles[i].ay -= source.m/p.m * prefac * dy;
        particles[i].az -= source.m/p.m * prefac * dz;
        
    }
}

void rebx_quad_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    for (int i=0; i<N; i++){
        const double* const mu_effcentral = rebx_get_param(sim->extras, particles[i].ap, "mu_effcentral");
        if (mu_effcentral != NULL){
            const double* const Rcentral = rebx_get_param(sim->extras, particles[i].ap, "Rcentral");
            if (Rcentral != NULL){
                rebx_calculate_quad_force(sim, particles, N, *mu_effcentral,*Rcentral, i); // only calculates force if a particle has mu_effcentral and Rcentral set
            }
        }
    }
}

static double rebx_calculate_quad_force_potential(struct reb_simulation* const sim, const double mu_eff, const double R, const int source_index){
    const struct reb_particle* const particles = sim->particles;
    const int _N_real = sim->N - sim->N_var;
    const struct reb_particle source = particles[source_index];
    double H = 0.;
    for (int i=0;i<_N_real;i++){
        if(i != 0){
            continue;
        }
        const struct reb_particle p = particles[i];
        const double dx = p.x - source.x;
        const double dy = p.y - source.y;
        const double dz = p.z - source.z;
        const double r2 = dx*dx + dy*dy + dz*dz;
        
        const double gamma = -4.0;
        const double G = sim->G;
        const double A = 3 * G * p.m / (4 * source.m);

        H -= p.m*mu_eff*A*pow(R,2)*pow(r2, (gamma+1.)/2.)/(gamma+1.);
        
    }
    return H;
}

double rebx_quad_force_potential(struct rebx_extras* const rebx){
    if (rebx->sim == NULL){
        rebx_error(rebx, ""); // rebx_error gives meaningful err
        return 0;
    }
    struct reb_simulation* sim = rebx->sim;
    const int N_real = sim->N - sim->N_var;
    struct reb_particle* const particles = sim->particles;
    double Htot = 0.;
    for (int i=0; i<N_real; i++){
        const double* const mu_effcentral = rebx_get_param(rebx, particles[i].ap, "mu_effcentral");
        if (mu_effcentral != NULL)
        {
            const double* const Rcentral = rebx_get_param(rebx, particles[i].ap, "Rcentral");
            if (Rcentral != NULL){
                Htot += rebx_calculate_quad_force_potential(sim, *mu_effcentral, *Rcentral, i);
            }
        }
    }
    return Htot;
}

// Add changing Sun mass requires dynamic A. We can use something like.

//double rebx_central_force_Acentral(const struct reb_particle p, const struct reb_particle primary, const double pomegadot, const double gamma){
//    struct reb_simulation* sim = p.sim;
//    const double G = sim->G;
//    const struct reb_orbit o = reb_orbit_from_particle(G, p, primary);
//    if (fabs(gamma+2.) < DBL_EPSILON){  // precession goes to 0 at r^-2, so A diverges for gamma=-2
//        reb_simulation_error(sim, "Precession vanishes for force law varying as r^-2, so can't initialize Acentral from a precession rate for gamma=-2)\n");
//        return 0.;
//    }
//    return G*primary.m*pomegadot/(1.+gamma/2.)/pow(o.d, gamma+2.)/o.n;
//}


