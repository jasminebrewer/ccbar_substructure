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

  int label = stoi(argv[2]); 
  
  // Generator. Process selection. 
  Pythia pythia;
  globalAnalysis analysis(pythia, static_cast<std::string>(argv[1]));
  analysis.initialize_pythia(label);

  vector<std::string> histogram_names {"eec2_all_"+to_string(label), "eec2_part_"+to_string(label), "eec2_part_cc_"+to_string(label), "eec2_part_bb_"+to_string(label), "eec3_all_"+to_string(label), "eec3_part_"+to_string(label), "eec2_med_all_"+to_string(label), "eec2_med_part_"+to_string(label), "eec3_med_all_"+to_string(label), "eec3_med_part_"+to_string(label)};
  analysis.declare_histograms(histogram_names);

  Histogram tmp(analysis._track_cuts.eecMin, analysis._track_cuts.eecMax, analysis._track_cuts.eec_bin_size);

  // set global values for splitting parameters
  struct splitting_params params(analysis._medium_params);

  // set up FastJet jet finder, with specified clustering algorithm, both for jet finding and reclustering 
  JetDefinition jet_def(analysis._jet_algo, analysis._track_cuts.jetR);
  JetDefinition jet_def_recl(analysis._jet_recl_algo, 5.*analysis._track_cuts.jetR, fastjet::E_scheme, fastjet::Best); // use larger radius for the reclustering to make sure all jet constituents are included, even if the jet area is irregular

  bool found_splitting;
  double med_weight=0.;

  // track total number of jets, and total number of ccbar-tagged jets
  int Njets=0;
  int NHFjets=0;

  bool jets_lose_energy = true;
  PseudoJet jet;
  
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
        
    for (int i=0; i<evt._jets.size(); i++) {

      Njets++;

      jet = evt._jets[i];
      double unmod_jet_pt = 0.0;
      if (evt._unmodified_jets.size()==evt._jets.size()) unmod_jet_pt = evt._unmodified_jets[i].perp();

      if (analysis._is_inclusive) {
        
        analysis._error_log << jet.perp() << ", "<< unmod_jet_pt << endl;
        tmp.set_to_zero();
	      EEC2(evt._tagged_cquarks, evt._tagged_cquarks, jet.pt(), tmp);
	      analysis._histograms.at("eec2_part_cc_"+to_string(label)).add(tmp._values);

	      tmp.set_to_zero();
	      EEC2(evt._tagged_bquarks, evt._tagged_bquarks, jet.pt(), tmp);
	      analysis._histograms.at("eec2_part_bb_"+to_string(label)).add(tmp._values);
      }

      if (!analysis._is_inclusive) {

  evt._has_pair = evt.get_pair(jet);
	if (!evt._has_pair) continue; // continue if the event doesn't have a particle/anti-particle pair

	//Splitting hardest_split = evt.find_hardest_splitting(jet);
	//cout << "hardest split: " << evt._py_event[hardest_split._in_index].id() << " -> " << evt._py_event[hardest_split._out1_index].id() << ", " << evt._py_event[hardest_split._out2_index].id() << endl;

	// evt.find_splitting();
	found_splitting = evt.find_splitting_v2(jet, analysis._track_cuts.jetR);
	if (!found_splitting) continue;

  if (analysis._match_splitting) {
	  // when using find_splitting_v2, can choose to overwrite the maxpt tagged particles with the final-state version of the particles found by the function
	  evt._maxpt_tagged_particle = evt.follow_to_final_state(evt._splitting._out1);
	  evt._maxpt_tagged_antiparticle = evt.follow_to_final_state(evt._splitting._out2);

	  // make sure the jet contains these particles as well
	  if (!contains(jet, evt._maxpt_tagged_particle) || !contains(jet, evt._maxpt_tagged_antiparticle)) continue;
  }
  
  NHFjets++;

	// use the parameters from the pythia event to compute the weight of the event for the medium modification
	params.z = evt._splitting._z;
	params.Eg = evt._splitting._Eg;
	params.pt2 = pow(evt._splitting._kt, 2.);
	med_weight = compute_medium_weight(&params, true); // gauss integration
            
	tmp.set_to_zero();
	EEC2(evt._tagged_particles, evt._tagged_particles, jet.pt(), tmp);
	
	analysis._error_log << jet.perp() << ", " << unmod_jet_pt << ", " << evt._maxpt_tagged_particle.perp() << ", " << evt._maxpt_tagged_antiparticle.perp() << ", " << (evt._maxpt_tagged_particle + evt._maxpt_tagged_antiparticle).perp() << ", " << med_weight;
  for (auto val: tmp._values) analysis._error_log << ", " << val;
  analysis._error_log << endl;

  if (med_weight<0. || med_weight>5.) continue;

  analysis._histograms.at("eec2_part_"+to_string(label)).add(tmp._values);
	analysis._histograms.at("eec2_med_part_"+to_string(label)).add(tmp._values, med_weight);

	tmp.set_to_zero();
	EEC3(evt._tagged_particles, evt._tagged_particles, jet.constituents(), jet.pt(), tmp);
	analysis._histograms.at("eec3_part_"+to_string(label)).add(tmp._values);
	analysis._histograms.at("eec3_med_part_"+to_string(label)).add(tmp._values, med_weight); 
      }// FOR FUTURE ANALYSES, SHOULD REALLY PUT "ALL" CORRELATORS INSIDE THIS LOOP

      // generic EEC: loops over all constituents, with the jet pt as the hard scale
      tmp.set_to_zero();
      EEC2(jet.constituents(), jet.constituents(), jet.pt(), tmp);
      analysis._histograms.at("eec2_all_"+to_string(label)).add(tmp._values);
      if (!analysis._is_inclusive) analysis._histograms.at("eec2_med_all_"+to_string(label)).add(tmp._values, med_weight);

      tmp.set_to_zero();
      EEC3(jet.constituents(), jet.constituents(), jet.constituents(), jet.pt(), tmp);
      analysis._histograms.at("eec3_all_"+to_string(label)).add(tmp._values);
      if (!analysis._is_inclusive) analysis._histograms.at("eec3_med_all_"+to_string(label)).add(tmp._values, med_weight);

    }// end of jet loop

  }// end of event loop. 
  // statistics
  analysis._pythia.stat();

  //analysis.normalize_histograms();

  // for (auto& histogram_pair : analysis._histograms) {
  //   histogram_pair.second.multiply(1.0/Njets);
  // }

  analysis.write_histograms();

  analysis._error_log << "number of jets: " << Njets << ", number of HF jets: " << NHFjets << endl;

  return 0;
}

