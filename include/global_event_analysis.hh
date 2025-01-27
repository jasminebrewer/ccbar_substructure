#ifndef __GLOBAL_EVENT_ANALYSIS_HH__
#define __GLOBAL_EVENT_ANALYSIS_HH__

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "Pythia8/Pythia.h"
#include "histograms.hh"
#include "constants.hh"
#include <map>

using namespace Pythia8;
using namespace fastjet;
using namespace std;


// object containing track, jet, and heavy flavor cuts and acceptances for all events
struct trackCuts {
  double trackEtaCut = 2.0;
  double trackPtMin = 0.0;
  double jetR = 0.5;
  double HFPtMin = 0.0;
  double JetPtMin = 100.0;
  double JetPtMax = 1e10;
  double JetEtaMin = -2.0;
  double JetEtaMax = 2.0;
  double eecMin = -2.5;
  double eecMax = -0.25;
  double eec_bin_size = 0.055;
  // double eecMin = 0.0;
  // double eecMax = 0.5;
  // double eec_bin_size = 0.02;
  double zcut=0.0;
  double beta=0.0;
};

struct mediumParams {
  double qL = 4.0;
  double L = 4.0 / constants::invGeVtofm;
  double omega_c = 60.; // GeV
  double n = 6.; // power law for the spectrum
  double T = 0.3; // temperature; GeV
  double alpha_med = 0.1; // medium coupling
  double mc2; // squared quark mass
  double jetR;
  string resolutionMode;
};


// class to extend the info in pseudoJet to include particle id and global event index information
class ExtraInfo : public PseudoJet::UserInfoBase {
public:
  // construct the user info f
  ExtraInfo(int pdg_id, int global_index=-1) {
    _pdg_id = pdg_id;
    _global_index = global_index;
  }
  // returns the pdg id of the particle
  int pdg_id() const {return _pdg_id;}

  // returns the global index of the particle in the event
  int global_index() const {return _global_index;}

private:
  int _pdg_id;
  int _global_index;
};



class globalAnalysis {
public:
  // copy constructor
  globalAnalysis(globalAnalysis& analysis) : _pythia(analysis._pythia), _jet_algo(analysis._jet_algo), _jet_recl_algo(analysis._jet_recl_algo), _track_cuts(analysis._track_cuts), _is_parton_level(analysis._is_parton_level), _is_inclusive(analysis._is_inclusive),
					     _particle_ids(analysis._particle_ids), _parton_ids(analysis._parton_ids), _medium_params(analysis._medium_params), _do_energy_loss(analysis._do_energy_loss), _recursive_daughters(analysis._recursive_daughters) {}

  // constructor a globalAnalysis instance from a pythia event and associated parameter file
  globalAnalysis(Pythia& pythia, string paramfile_name) : _pythia(pythia) {
    _parameter_file.open(paramfile_name.c_str(), ios::in);
    assert(_parameter_file.is_open() && "cannot open the specified parameter file.");
  }

  // functions
  void initialize_pythia(int label);
  void declare_histograms(vector<string> histogram_names);
  void normalize_histograms();
  void write_histograms();

  // members
  Pythia& _pythia;                      // pythia instance
  int _n_events;                        // number of events to generate

  string _file_label;                   // label to append to the histogram file unique for each run (e.g., inc, qq, cc, or bb)
  ofstream _error_log;                  // log file for errors
  fstream _parameter_file;              // file from which to read global event analysis parameters

  JetAlgorithm _jet_algo;               // algorithm for the jet clustering
  JetAlgorithm _jet_recl_algo;          // algorithm for the reclustering step
  string _FC_mode;                      // string specifying the type of flavor cone method to use
  bool _match_splitting;                // boolean specifying whether or not to match particles between pythia splitting and in the reclustering
  bool _apply_sd_to_all;                // boolean for whether to apply the softdrop cut to the whole jet, or just for the purpose of event selection
  bool _recursive_daughters;            // boolean for whether to determine the "MC truth" based on just the splitting duaghters themselves (false) or recursively including all subsequent daughters (true)

  trackCuts _track_cuts;                // struct containing kinematic cuts for the analysis
  mediumParams _medium_params;          // struct containing parameters for the medium modification

  bool _do_energy_loss;                 // flag for whether to compute jet energy loss
  bool _is_parton_level;                // flag for whether the event is parton- or hadron-level
  string _hadron_mode;                  // flag for whether to identify only D0/D0bar ('D/B') or all hadrons (or equivalently for B) ('all')
  bool _is_inclusive;                   // whether to tag particular particles in the event or consider all equally. If false, particle_ids should not be empty
  vector<int> _particle_ids;            // particle ids of particles to tag (e.g., c and cbar)
  vector<int> _parton_ids;              // for the case that hadronization is on, also set the ids of associated partons

  map<string, Histogram> _histograms;   // histogram object and a string handle used internally to identify each histogram

};

#endif // __GLOBAL_EVENT_ANALYSIS_HH__
