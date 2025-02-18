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

  // set up FastJet jet finder, with specified clustering algorithm, both for jet finding and reclustering 
  JetDefinition jet_def(analysis._jet_algo, analysis._track_cuts.jetR);

  bool jets_lose_energy = true;
  PseudoJet jet;
  
  // Begin event loop
  for (int iEvent = 0; iEvent < analysis._n_events; ++iEvent) {

    if (!analysis._pythia.next()) continue;

    EventCCbar evt(analysis);
    evt._jet_def = jet_def;
    analysis._is_inclusive = true;

    // read in particles from the pythia event
    evt.read_event(); 

    // if inclusive, consider up to two jets that pass the cuts in an event
    // if not inclusive, consider up to two jets that pass the cuts and consider the one that has at least one pair of particles with the indices of tagged particles
    evt.cluster_jets();

    if (evt._jets.size()==0) continue; // event has no jets passing the selections
        
    for (int i=0; i<evt._jets.size(); i++) {

      jet = evt._jets[i];
      double unmod_jet_pt = 0.0;
      if (evt._unmodified_jets.size()>=evt._jets.size()) unmod_jet_pt = evt._unmodified_jets[i].perp();

      analysis._error_log << jet.perp() << ", "<< unmod_jet_pt << endl;

    }// end of jet loop

  }// end of event loop. 
  // statistics
  analysis._pythia.stat();

  //analysis.normalize_histograms();

  // for (auto& histogram_pair : analysis._histograms) {
  //   histogram_pair.second.multiply(1.0/Njets);
  // }

  // analysis.write_histograms();

  // analysis._error_log << "number of jets: " << Njets << ", number of HF jets: " << NHFjets << endl;
  // cout << "number of jets: " << Njets << ", number of HF jets: " << NHFjets << ", avg pt HF: " << avgptHF << endl;

  return 0;
}

