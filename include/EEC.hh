#ifndef __EEC_HH__
#define __EEC_HH__

#include "fastjet/PseudoJet.hh"
#include "histograms.hh"

using namespace fastjet;

// function to compute the 2-point energy correlator between particles1 and particles2
void EEC2(vector<PseudoJet> particles1, vector<PseudoJet> particles2, double pT_scale, Histogram& hist);

// function to compute the 3-point energy correlator between particles1, particles2, and particles3
void EEC3(vector<PseudoJet> particles1, vector<PseudoJet> particles2, vector<PseudoJet> particles3, double pT_scale, Histogram& hist);

#endif // __EEC_HH__
