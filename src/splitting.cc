#include "fastjet/PseudoJet.hh"
#include "splitting.hh"
#include "global_event_analysis.hh"


using namespace fastjet;


/**
 * @brief function to compute the (3-vector) dot product between the momenta of two jets
 * 
 * @param jet1, jet2: PseudoJet objects whose momenta you want to dot product
 * @return 3-vector momentum dot product
*/
double dot(PseudoJet jet1, PseudoJet jet2) {

  return jet1.px()*jet2.px() + jet1.py()*jet2.py() + jet1.pz()*jet2.pz();

}


/**
 * @brief function to get the momentum sharing fraction z between two subjets that are the result of a splitting
 * 
 * @param leading_subjet: higher-pt subjet of a splitting
 * @param subleading_subjet: lower-pt subjet of a splitting
 * @return momentum fraction carried by the subleading subjet
*/
double get_z( PseudoJet leading_subjet, PseudoJet subleading_subjet) {

  PseudoJet sum = leading_subjet + subleading_subjet;

  return dot(sum, subleading_subjet) / dot(sum, sum);
}

/**
 * @brief function to get the tranvserse momentum associated with a splitting
 * 
 * @param leading_subjet: higher-pt subjet of a splitting
 * @param subleading_subjet: lower-pt subjet of a splitting
 * @return transverse momentum of the splitting
*/double get_kt( PseudoJet leading_subjet, PseudoJet subleading_subjet) {

  PseudoJet sum = leading_subjet + subleading_subjet;

  double dot_prod, norm2_jet, norm2_subjet, proj_norm;
  norm2_jet = dot(sum, sum);
  norm2_subjet = dot(leading_subjet, leading_subjet);
  dot_prod = dot(sum, leading_subjet);
  proj_norm = norm2_subjet - (dot_prod*dot_prod)/norm2_jet;
  
  return sqrt(proj_norm);
}

/**
 * @brief function to set the physical values of a splitting (Eg, z, etc) given that the PseudoJets associated with the splitting
 * itself (in, out1, and out2) have been defined.
*/
void Splitting::set_values() {

  _Eg = _in.e();
  // _Egv2 = (_out1+_out2).e();
  _kt = get_kt(_out1, _out2);
  _z = get_z(_out1, _out2);
  _dR = _out1.delta_R(_out2);
  _pt = _in.pt();
  _virt = sqrt( (pow(_kt,2.0) + pow(constants::CHARM_MASS,2.0)) / (_z*(1.0-_z)) );
  _virt_v2 = (_out1+_out2).m();
}
