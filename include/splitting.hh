#ifndef __SPLITTING_HH__
#define __SPLITTING_HH__

#include "fastjet/PseudoJet.hh"

using namespace fastjet;

class Splitting {
public:

  // default constructor: fill the splitting with empty PseudoJets
  Splitting() : _in(PseudoJet()), _out1(PseudoJet()), _out2(PseudoJet()) {}
  
  Splitting(PseudoJet in, PseudoJet out1, PseudoJet out2) :
    _in(in), _out1(out1), _out2(out2) {}
  
  // take a splitting where the _in, _out1, and _out2 are defined and calculate the associated kinematics of the splitting (Eg, z, etc)
  void set_values();

  fastjet::PseudoJet _in;    // mother ("in") particle of the splitting
  fastjet::PseudoJet _out1;  // higher-pt daughter particle of the splitting
  fastjet::PseudoJet _out2;  // lower-pt daughter particle of the splitting
  bool _is_valid;            // false if the splitting is not valid; e.g. if a splitting is not found in pythia
  int _level;                // integer for the level in the shower the splitting occured
  bool _is_primary;          // flag for the whether the splitting is on the primary Lund plane, i.e. always on the higher-pt branch

  double _Eg;                // energy of the mother gluon
  double _pt;                // pt of the mother gluon
  double _kt;                // kt, z, dR of the splitting
  double _z;
  double _dR;
  double _virt;              // virtuality of the splitting defined from kt and z
  double _virt_v2;           // virtuality of the splitting defined from subjet invariant mass
};

#endif // __SPLITTING_HH__