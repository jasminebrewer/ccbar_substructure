#ifndef __MEDIUM_MOD_HH__
#define __MEDIUM_MOD_HH__

struct splitting_params {
    double z; 
    double Eg; 
    double mc2; 
    double pt2; 
    double qL; 
    double L;
    
    splitting_params(mediumParams params) : mc2(params.mc2), qL(params.qL), L(params.L) {}
}; 

double compute_medium_weight(void * p, bool gauss_integrator);

#endif // __MEDIUM_MOD_HH__
