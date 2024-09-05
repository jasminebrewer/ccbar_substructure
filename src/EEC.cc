#include "fastjet/PseudoJet.hh"
#include "histograms.hh"

using namespace fastjet;

/**
 * @brief function to compute the two-point energy correlator between particles1 and particles2, and store the result in a histogram
 * 
 * @param particles1, particles2: vector of PseudoJet objects (particles) to take the two-point correlator of
 * @param pT_scale: pT scale entering the energy correlator (generally, the jet energy)
 * @param hist: Histogram object to store the resulting correlator in
 */
void EEC2(vector<PseudoJet> particles1, vector<PseudoJet> particles2, double pT_scale, Histogram& hist) {

  double x_2L, z1, z2;
  for (auto p1: particles1) {
    z1 = p1.pt() / pT_scale;
    for (auto p2: particles2) {
      z2 = p2.pt() / pT_scale;
      x_2L = p1.delta_R(p2);

      hist.fill(log10(x_2L), z1*z2);
    }
  }
}

/**
 * @brief function to compute the three-point energy correlator between particles1, particles2, and particles3, and store the result in a histogram
 * 
 * @param particles1, particles2, particles3: vector of PseudoJet objects (particles) to take the two-point correlator of
 * @param pT_scale: pT scale entering the energy correlator (generally, the jet energy)
 * @param hist: Histogram object to store the resulting correlator in
 */
void EEC3(vector<PseudoJet> particles1, vector<PseudoJet> particles2, vector<PseudoJet> particles3, double pT_scale, Histogram& hist) {

  double x_2L, z1, z2, z3;
  for (auto p1: particles1) {
    z1 = p1.pt() / pT_scale;
    for (auto p2: particles2) {
      z2 = p2.pt() / pT_scale;
      for (auto p3: particles3) {
        z3 = p3.pt() / pT_scale;
        x_2L = p1.delta_R(p2);

	      hist.fill(log10(x_2L), z1*z2*z3 );
      }
    }
  }
}
