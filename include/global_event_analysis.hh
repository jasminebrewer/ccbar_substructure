#ifndef __GLOBAL_EVENT_ANALYSIS_HH__
#define __GLOBAL_EVENT_ANALYSIS_HH__

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "Pythia8/Pythia.h"
#include "histograms.hh"
#include <map>

using namespace Pythia8;
using namespace fastjet;
using namespace std;

namespace constants
{
  // identifying particle ids used by pythia
  const int CHARM = 4;
  const int ANTICHARM = -4;
  const int BOTTOM = 5;
  const int ANTIBOTTOM = -5;
  const int DOWN = 1;
  const int ANTIDOWN = -1;
  const int GLUON = 21;
  const int D0 = 421;
  const int D0BAR = -421;
  // particles masses
  const double CHARM_MASS = 1.27;
  const double BOTTOM_MASS = 4.18;
  // status code for particles in pythia
  const int INCOMING_HARD = -21; // incoming from hard process
  const int INCOMING_SUB = -31; // incoming from subprocess
  const int INCOMING_ISR = -41; // incoming from initial state radiation
  const int OUTGOING_HARD = -23; // outgoing from hard process
  const int OUTGOING_SUB = -33; // outgoing from subprocess
  const int OUTGOING_ISR = -43; // outgoing from initial state radiation
  const int SHOWER = -51; // produced by the shower
  // unit conversions
  const double invGeVtofm = 0.1973; // conversion for length units in GeV^{-1} to fermi
}


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
					     _particle_ids(analysis._particle_ids), _qhatL(analysis._qhatL), _L(analysis._L), _mc2(analysis._mc2) {}

  // constructor a globalAnalysis instance from a pythia event and associated parameter file
  globalAnalysis(Pythia& pythia, string paramfile_name) : _pythia(pythia) {
    _parameter_file.open(paramfile_name.c_str(), ios::in);
    assert(_parameter_file.is_open() && "cannot open the specified parameter file.");
  }

  // functions
  void initialize_pythia();
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

  trackCuts _track_cuts;                // struct containing kinematic cuts for the analysis

  bool _is_parton_level;                // flag for whether the event is parton- or hadron-level
  bool _is_inclusive;                   // whether to tag particular particles in the event or consider all equally. If false, particle_ids should not be empty
  vector<int> _particle_ids;            // particle ids of particles to tag (e.g., c and cbar)

  map<string, Histogram> _histograms;   // histogram object and a string handle used internally to identify each histogram

  // global parameters for splittings and their medium modification
  double _qhatL;        // strength of the medium interaction (qhat) times length
  double _L;            // length of the medium
  double _mc2;          // quark mass

};

#endif // __GLOBAL_EVENT_ANALYSIS_HH__
