/**
 * Post-Newtonian correction from general relativity
 * 
 * This example shows how to add post-newtonian corrections to REBOUND simulations with reboundx.
 * If you have GLUT installed for the visualization, press 'w' and/or 'c' for a clearer view of
 * the whole orbit.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]){
	struct reb_simulation* sim = reb_create_simulation();
	// Setup constants
	sim->dt 		= 1.e-2;		// timestep.
	sim->integrator	= REB_INTEGRATOR_WHFAST;
	//sim->integrator	= REB_INTEGRATOR_IAS15;

	struct reb_particle p = {0}; 
	p.m  	= 1.;	
	reb_add(sim, p); 

	double m = 0.;
	double a = 1.; // put planet close to enhance precession (this would put planet inside the Sun!)
	double e = 0.2;
	double omega = 0.;
	double f = 0.;

	struct reb_particle p1 = reb_tools_orbit2d_to_particle(sim->G, p, m, a, e, omega, f);
	reb_add(sim,p1);
	reb_move_to_com(sim);
	
	rebx_init(sim); // initialize reboundx
	double c = C_DEFAULT; // Have to set the speed of light in appropriate units (set by G and your initial conditions).  Here we use the value in default units of AU/(yr/2pi)	
	rebx_add_gr(sim,c); // add postnewtonian correction.  


	double tmax = 1.e2;
	reb_integrate(sim, tmax); 
}
