#include "fastjet/PseudoJet.hh"

using namespace fastjet;

PseudoJet scaleMomentum(PseudoJet particle, double scale);
vector<PseudoJet> compute_jet_modification( vector<PseudoJet> jet_constituents, void * p );