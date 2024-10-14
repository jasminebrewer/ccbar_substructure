#ifndef __CCBAR_ANALYSIS_HH__
#define __CCBAR_ANALYSIS_HH__

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "Pythia8/Pythia.h"
#include "histograms.hh"
#include "splitting.hh"
#include <map> 

#include "global_event_analysis.hh"

using namespace Pythia8;
using namespace fastjet;
using namespace std;


//=========================================================================
/// @brief class that holds information about events needed for ccbar-tagged analysis
class EventCCbar : globalAnalysis {
public:

  /// @brief constructor for ccbar events
  EventCCbar(globalAnalysis& analysis) : globalAnalysis(analysis) , _py_event(analysis._pythia.event) {}

  void read_event();
  void cluster_jets();
  void find_splitting();
  bool find_splitting_v2(PseudoJet jet, double jetR);
  void calculate_splitting_level();
  PseudoJet follow_to_final_state(PseudoJet particle);
  Splitting do_iterative_reclustering(PseudoJet jet);
  Splitting do_flavor_cone(string FC_mode);
  Splitting get_random_splitting(PseudoJet jet);
  bool get_pair();
  bool get_pair(PseudoJet jet);

  //Splitting find_hardest_splitting(PseudoJet jet);

  // members

  Event& _py_event;                     // the input pythia event
  vector<PseudoJet> _final_particles;   // vector of all final state particles in the event

  JetDefinition _jet_def;               // the jet definition
  ClusterSequence _cluster_seq;         // cluster sequence from fastjet (necessary to access jet substructure)

  JetDefinition _jet_def_recl;          // jet definition for the reclustering step

  vector<PseudoJet> _jets;              // up to two jets in the event passing the cuts
  vector<PseudoJet> _unmodified_jets;

  vector<PseudoJet> _tagged_particles;  // vector of PseudoJets containing the tagged final-state particles
  // temporary addition!!
  vector<PseudoJet> _tagged_cquarks;
  vector<PseudoJet> _tagged_bquarks;

  vector<PseudoJet> _tagged_partons;    // for the case where hadronization is on, also store the tagged particles at parton-level
  PseudoJet _maxpt_tagged_antiparticle; // PseudoJet containing the maximum pt tagged particle in the event
  PseudoJet _maxpt_tagged_particle;     // PsuedoJet containing the maximum pt antiparticle in the event
  bool _has_pair;                       // true if the event contains a particle/ antiparticle pair of the tagged type; false otherwise

  Splitting _splitting;                 // the MC-level splitting associated with the tagged particles
  Splitting _recl_splitting;            // the splitting, according to iterative reclustering, between the tagged particles
};


#endif // __CCBAR_ANALYSIS_HH__
