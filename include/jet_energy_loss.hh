#include "fastjet/PseudoJet.hh"

using namespace fastjet;

struct energy_loss_params {double qL; double L; double omega_c; double n; double T; double alpha_med; double jetR;}; 

double estimate_energy_loss( PseudoJet jet, void * p );