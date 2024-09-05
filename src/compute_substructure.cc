#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "ccbar_analysis.hh"
#include "complex_Ei.hh"
#include "histograms.hh"
#include "medium_mod.hh"
#include <time.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <complex>

using namespace Pythia8;
using namespace fastjet;
using namespace std;

bool contains(fastjet::PseudoJet& jet, fastjet::PseudoJet particle);

int main(int argc, char* argv[]) {

  if (argc!=3) {
    cout << "Incorrect number of command line arguments -- command line argument should be the name of a parameter file" << endl;
    return -1;
  }

  // Generator. Process selection. 
  Pythia pythia;
  globalAnalysis analysis(pythia, static_cast<std::string>(argv[1]));
  analysis.initialize_pythia(stoi(argv[2]));

  // set global values for parameters needed to compute the medium modification of the splitting
  struct splitting_params params(analysis._medium_params);

  // set up FastJet jet finder, with specified clustering algorithm, both for jet finding and reclustering 
  JetDefinition jet_def(analysis._jet_algo, analysis._track_cuts.jetR);
  JetDefinition jet_def_recl(analysis._jet_recl_algo, 5.*analysis._track_cuts.jetR, fastjet::E_scheme, fastjet::Best); // use larger radius for the reclustering to make sure all jet constituents are included, even if the jet area is irregular


  double z_cut = analysis._track_cuts.zcut;
  double beta  = analysis._track_cuts.beta;
  fastjet::contrib::SoftDrop sd(beta, z_cut);
  cout << "SoftDrop groomer is: " << sd.description() << endl;

  bool found_splitting;
  double med_weight=0.;

  Splitting recl_splitting, mod_recl_splitting, splitting_cc;
  
  // Begin event loop
  for (int iEvent = 0; iEvent < analysis._n_events; ++iEvent) {

    if (!analysis._pythia.next()) continue;

    EventCCbar evt(analysis);
    evt._jet_def = jet_def;
    evt._jet_def_recl = jet_def_recl;

    // read in particles from the pythia event
    evt.read_event();

    // if inclusive, consider up to two jets that pass the cuts in an event
    // if not inclusive, consider up to two jets that pass the cuts and consider the one that has at least one pair of particles with the indices of tagged particles
    // if do_energy_loss is true, first modify the particles and then perform jet selections
    evt.cluster_jets();

    if (evt._jets.size()==0) continue; // event has no jets passing the selections
        
    for (auto jet: evt._jets) {

      if (!analysis._is_inclusive) {

	if (!evt._has_pair) continue; // continue if the event doesn't have a particle/anti-particle pair


	// do soft-drop grooming on the jet to remove soft splittings                                                                                                  
	fastjet::PseudoJet sd_jet = sd(jet);

	// reject jets that do not contain both the maxpt particle and the maxpt antiparticle after the softdrop grooming
	if (! (contains(sd_jet, evt._maxpt_tagged_particle) && contains(sd_jet, evt._maxpt_tagged_antiparticle)) ) continue;

	//Splitting hardest_split = evt.find_hardest_splitting(jet);
	//cout << "hardest split: " << evt._py_event[hardest_split._in_index].id() << " -> " << evt._py_event[hardest_split._out1_index].id() << ", " << evt._py_event[hardest_split._out2_index].id() << endl;

	// evt.find_splitting();
  cout<< "finding splitting" << endl;
	found_splitting = evt.find_splitting_v2(jet, analysis._track_cuts.jetR);
	if (!found_splitting) continue;

	// when using find_splitting_v2, overwrite the maxpt tagged particles with the final-state version of the particles found by the function
	evt._maxpt_tagged_particle = evt.follow_to_final_state(evt._splitting._out1);
	evt._maxpt_tagged_antiparticle = evt.follow_to_final_state(evt._splitting._out2);

	// make sure the jet contains these particles as well
	if (!contains(jet, evt._maxpt_tagged_particle) || !contains(jet, evt._maxpt_tagged_antiparticle)) continue;

cout << "level" << endl;
	evt.calculate_splitting_level(); // calculate the level of the splitting found above

cout << "iterative reclustering" << endl;
	recl_splitting = evt.do_iterative_reclustering(jet); // recluster the jet

// cout << "flavor cone" << endl;
//   splitting_cc = evt.do_flavor_cone(analysis._FC_mode);


	// use the parameters from the pythia event to compute the weight of the event for the medium modification
	params.z = evt._splitting._z;
	params.Eg = evt._splitting._Eg;
	params.pt2 = pow(evt._splitting._kt, 2.);
  cout << "medium weight" << endl;
	med_weight = compute_medium_weight(&params, true); // gauss integration

analysis._error_log << med_weight << ", " << evt._splitting._is_valid << ", " << evt._splitting._level << ", " << evt._splitting._Eg << ", " << evt._splitting._pt << ", " << evt._splitting._kt << ", " << evt._splitting._z << ", " << evt._splitting._dR << ", " << evt._splitting._virt << ", " << recl_splitting._is_primary << ", " << recl_splitting._level << ", " << recl_splitting._Eg << ", " << recl_splitting._pt << ", " << recl_splitting._kt << ", " << recl_splitting._z << ", " << recl_splitting._dR << ", " << recl_splitting._virt << endl;
  // analysis._error_log << med_weight << ", " << evt._splitting._is_valid << ", " << evt._splitting._level << ", " << evt._splitting._Eg << ", " << evt._splitting._pt << ", " << evt._splitting._kt << ", " << evt._splitting._z << ", " << evt._splitting._dR << ", " << evt._splitting._virt << ", " << recl_splitting._is_primary << ", " << recl_splitting._level << ", " << recl_splitting._Eg << ", " << recl_splitting._pt << ", " << recl_splitting._kt << ", " << recl_splitting._z << ", " << recl_splitting._dR << ", " << recl_splitting._virt << ", " << splitting_cc._Eg << ", " << splitting_cc._pt << ", " << splitting_cc._kt << ", " << splitting_cc._z << ", " << splitting_cc._dR << ", " << splitting_cc._virt << endl;
  // analysis._error_log << jet.perp() << ", " << modified_jet.perp() << ", " << med_weight << ", " << evt._splitting._is_valid << ", " << evt._splitting._level << ", " << evt._splitting._Eg << ", " << evt._splitting._pt << ", " << evt._splitting._kt << ", " << evt._splitting._z << ", " << evt._splitting._dR << ", " << evt._splitting._virt << ", " << recl_splitting._is_primary << ", " << recl_splitting._level << ", " << recl_splitting._Eg << ", " << recl_splitting._pt << ", " << recl_splitting._kt << ", " << recl_splitting._z << ", " << recl_splitting._dR << ", " << recl_splitting._virt << ", " << mod_recl_splitting._level << ", " << mod_recl_splitting._Eg << ", " << mod_recl_splitting._pt << ", " << mod_recl_splitting._kt << ", " << mod_recl_splitting._z << ", " << mod_recl_splitting._dR << ", " << mod_recl_splitting._virt << ", " << splitting_cc._Eg << ", " << splitting_cc._pt << ", " << splitting_cc._kt << ", " << splitting_cc._z << ", " << splitting_cc._dR << ", " << splitting_cc._virt << endl;
      }

    }// end of jet loop

  }// end of event loop. 
  // statistics
  analysis._pythia.stat();

  // analysis.normalize_histograms();


  // analysis.write_histograms();

  return 0;
}

