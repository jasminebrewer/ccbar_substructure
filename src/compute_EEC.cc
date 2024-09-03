#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "ccbar_analysis.hh"
#include "EEC.hh"
#include "complex_Ei.hh"
#include "histograms.hh"
#include "medium_mod.hh"
#include "jet_energy_loss.hh"
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

  vector<std::string> histogram_names {"eec2_all", "eec2_part", "eec3_all", "eec3_part", "eec2_med_all", "eec2_med_part", "eec3_med_all", "eec3_med_part"};
  analysis.declare_histograms(histogram_names);

  Histogram tmp(analysis._track_cuts.eecMin, analysis._track_cuts.eecMax, analysis._track_cuts.eec_bin_size);

  // set global values for splitting parameters
  struct splitting_params params;
  params.qL = analysis._qhatL;
  params.L = analysis._L;
  params.mc2 = analysis._mc2;

  // set up FastJet jet finder, with specified clustering algorithm, both for jet finding and reclustering 
  JetDefinition jet_def(analysis._jet_algo, analysis._track_cuts.jetR);
  JetDefinition jet_def_recl(analysis._jet_recl_algo, 5.*analysis._track_cuts.jetR, fastjet::E_scheme, fastjet::Best); // use larger radius for the reclustering to make sure all jet constituents are included, even if the jet area is irregular

  bool found_splitting;
  double med_weight=0.;

  struct energy_loss_params eloss_params;
  eloss_params.qL = analysis._qhatL;
  eloss_params.L = analysis._L;
  eloss_params.jetR = analysis._track_cuts.jetR;
  eloss_params.omega_c = 60.; // GeV
  eloss_params.n = 6.;
  eloss_params.T = 0.3; // GeV
  eloss_params.alpha_med = 0.1;

  // track total number of jets, and total number of ccbar-tagged jets
  int Njets=0;
  int NHFjets=0;

  bool jets_lose_energy = true;
  
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
    evt.cluster_jets();

    //cout << "num jets: " << evt._jets.size() << endl;

    if (evt._jets.size()==0) continue; // event has no jets passing the selections
        
    for (auto jet: evt._jets) {

      Njets++;

      if (jets_lose_energy) {
      
        vector<PseudoJet> modified_jet_particles = compute_jet_modification(jet, &eloss_params);

        JetDefinition mod_jet_def(antikt_algorithm, 5.0*analysis._track_cuts.jetR); // make a very large jet to make sure to cluster all of the particles
        ClusterSequence mod_cs(modified_jet_particles, mod_jet_def);
        vector<PseudoJet> mod_jets = sorted_by_pt(mod_cs.inclusive_jets(0.0));
        PseudoJet modified_jet = mod_jets[0];
      }

      //if (analysis._is_inclusive) analysis._error_log << jet.perp() << ", " << modified_jet.perp() << endl;


      if (!analysis._is_inclusive) {

	if (!evt._has_pair) continue; // continue if the event doesn't have a particle/anti-particle pair

	NHFjets++;

	//Splitting hardest_split = evt.find_hardest_splitting(jet);
	//cout << "hardest split: " << evt._py_event[hardest_split._in_index].id() << " -> " << evt._py_event[hardest_split._out1_index].id() << ", " << evt._py_event[hardest_split._out2_index].id() << endl;

	// evt.find_splitting();
	found_splitting = evt.find_splitting_v2(jet, analysis._track_cuts.jetR);
	if (!found_splitting) continue;

	// when using find_splitting_v2, overwrite the maxpt tagged particles with the final-state version of the particles found by the function
	evt._maxpt_tagged_particle = evt.follow_to_final_state(evt._splitting._out1);
	evt._maxpt_tagged_antiparticle = evt.follow_to_final_state(evt._splitting._out2);

	// make sure the jet contains these particles as well
	if (!contains(jet, evt._maxpt_tagged_particle) || !contains(jet, evt._maxpt_tagged_antiparticle)) continue;

	// use the parameters from the pythia event to compute the weight of the event for the medium modification
	params.z = evt._splitting._z;
	params.Eg = evt._splitting._Eg;
	params.pt2 = pow(evt._splitting._kt, 2.);
	med_weight = compute_medium_weight(&params, true); // gauss integration
            
	tmp.set_to_zero();
	EEC2(evt._tagged_particles, evt._tagged_particles, jet.pt(), tmp);
	analysis._histograms.at("eec2_part").add(tmp._values);
	analysis._histograms.at("eec2_med_part").add(tmp._values, med_weight);

	tmp.set_to_zero();
	EEC3(evt._tagged_particles, evt._tagged_particles, jet.constituents(), jet.pt(), tmp);
	analysis._histograms.at("eec3_part").add(tmp._values);
	analysis._histograms.at("eec3_med_part").add(tmp._values, med_weight); 
      }

      // generic EEC: loops over all constituents, with the jet pt as the hard scale
      tmp.set_to_zero();
      EEC2(jet.constituents(), jet.constituents(), jet.pt(), tmp);
      analysis._histograms.at("eec2_all").add(tmp._values);
      if (!analysis._is_inclusive) analysis._histograms.at("eec2_med_all").add(tmp._values, med_weight);

      tmp.set_to_zero();
      EEC3(jet.constituents(), jet.constituents(), jet.constituents(), jet.pt(), tmp);
      analysis._histograms.at("eec3_all").add(tmp._values);
      if (!analysis._is_inclusive) analysis._histograms.at("eec3_med_all").add(tmp._values, med_weight);

    }// end of jet loop

  }// end of event loop. 
  // statistics
  analysis._pythia.stat();

  analysis.normalize_histograms();

  analysis.write_histograms();

  analysis._error_log << "number of jets: " << Njets << ", number of HF jets: " << NHFjets << endl;

  return 0;
}

