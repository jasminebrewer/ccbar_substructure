#ifndef __MEDIUM_MOD_HH__
#define __MEDIUM_MOD_HH__

#include "global_event_analysis.hh"

struct splitting_params {
    double z = 0.; 
    double Eg = 0.; 
    double mc2 = 0.; 
    double pt2 = 0.; 
    double qL = 0.; 
    double L = 0.;
    
    splitting_params(mediumParams params) : mc2(params.mc2), qL(params.qL), L(params.L) {}
    splitting_params() {}
}; 

double compute_medium_weight(void * p, bool gauss_integrator);

#endif // __MEDIUM_MOD_HH__
