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

  Histogram spec_hist_vac(0.0, 500.0, 10.0, "spectrum_vac");
  Histogram spec_hist_med(0.0, 500.0, 10.0, "spectrum_med");

  bool found_splitting;
  double med_weight=0.;

  Splitting recl_splitting, mod_recl_splitting, splitting_cc;

  PseudoJet jet, unmod_jet;
  
  // Begin event loop
  for (int iEvent = 0; iEvent < analysis._n_events; ++iEvent) {

    if (!analysis._pythia.next()) continue;

    EventCCbar evt(analysis);
    evt._jet_def = jet_def;

    // read in particles from the pythia event
    evt.read_event();

    // if inclusive, consider up to two jets that pass the cuts in an event
    // if not inclusive, consider up to two jets that pass the cuts and consider the one that has at least one pair of particles with the indices of tagged particles
    // if do_energy_loss is true, first modify the particles and then perform jet selections
    evt.cluster_jets();

    if (evt._jets.size()==0) continue; // event has no jets passing the selections
        
    for (int i=0; i<evt._jets.size(); i++) {

      jet = evt._jets[i];
      double unmod_jet_pt = 0.0;
      if (evt._unmodified_jets.size()==evt._jets.size()) unmod_jet_pt = evt._unmodified_jets[i].perp();

      if (analysis._is_inclusive) {
        
        analysis._error_log << jet.perp() << ", "<< unmod_jet_pt << ", 1" << endl;
      }

      if (!analysis._is_inclusive) {

        // check if the jet contains at least two of the tagged_particles
        evt._has_pair = evt.get_pair(jet);
        if (!evt._has_pair) continue;

	      found_splitting = evt.find_splitting_v2(jet, analysis._track_cuts.jetR, analysis._recursive_daughters);
	      if (!found_splitting) continue;

        if (analysis._match_splitting) {
	        // when using find_splitting_v2, can choose to overwrite the maxpt tagged particles with the final-state version of the particles found by the function
	        evt._maxpt_tagged_particle = evt.follow_to_final_state(evt._splitting._out1);
	        evt._maxpt_tagged_antiparticle = evt.follow_to_final_state(evt._splitting._out2);

	        // make sure the jet contains these particles as well
	        if (!contains(jet, evt._maxpt_tagged_particle) || !contains(jet, evt._maxpt_tagged_antiparticle)) continue;
        }

	      // use the parameters from the pythia event to compute the weight of the event for the medium modification
	      params.z = evt._splitting._z;
	      params.Eg = evt._splitting._Eg;
	      params.pt2 = pow(evt._splitting._kt, 2.);
	      med_weight = compute_medium_weight(&params, true); // gauss integration

        spec_hist_vac.fill( jet.perp(), 1.0);
        spec_hist_med.fill( jet.perp(), med_weight);

        analysis._error_log << jet.perp() << ", "<< unmod_jet_pt << ", " << med_weight << endl;
      }

    }// end of jet loop

  }// end of event loop. 
  // statistics
  analysis._pythia.stat();

  spec_hist_vac.write();
  spec_hist_med.write();

  analysis.write_histograms();

  return 0;
}

