#include <fstream> 
#include "global_event_analysis.hh"
#include "constants.hh"
#include <boost/property_tree/ptree.hpp> // used for processing ini initialization files
#include <boost/property_tree/ini_parser.hpp>


/** 
 * @brief helper function to determine if a particle has charm or bottom quarks from its particle id
*/
bool has_charm(Particle p) {

  int id = p.idAbs();
  return ((id / 100) % 10 == 4 || (id / 1000) % 10 == 4); // check for 4 in 100s or 1000s place
}

bool has_bottom(Particle p) {
  int id = p.idAbs();
  return ((id / 100) % 10 == 5 || (id / 1000) % 10 == 5); // check for 4 in 100s or 1000s place
}

/**
 * @brief function to set parameters of the analysis and initialize pythia from values specified in the parameter file
 */
void globalAnalysis::initialize_pythia(int label) {

    boost::property_tree::ptree tree;
    boost::property_tree::ini_parser::read_ini(_parameter_file, tree);

    set<int> chadrons = {411,421,10411,10421,413,423,10413,10423,20413,20423,415,425,431,10431,433,10433,20433,435,4122,4222,4212,4112,4224,4214,4232,4132,4322,4312,4324,4314,4332,4334,4412,4422,4414,4424,4432,4434,4444,441,10441,100441,443,10443,20443,100443,30443,445};
    set<int> bhadrons = {511,521,513,523,515,525,531,533,535,541,543,545,10511,10521,10513,10523,20513,20523,10531,10533,20533,10541,10543,20543,5122,5112,5212,5222,5114,5214,5224,5132,5232,5312,5322,5314,5324,5332,5334,5142,5242,5412,5422,5414,5424,5342,5432,5434,5442,5444,5512,5522,5514,5524,5532,5534,5542,5544,5554,551,10551,553,10553,20553,100553,200553,555};

    //* load parameters for the pythia generation *//

    _n_events = tree.get<int>("pythia.maxneventsperjob", 10000);
    _pythia.readString("Main:numberOfEvents = "+to_string(_n_events));

    set<string> allowed_processes = {"all", "gg"};  
    string process = tree.get<string>("selection.process", "all");
    if (allowed_processes.find(process) == allowed_processes.end()) {
        throw invalid_argument("Invalid process specified: " + process);
    }
    if (process == "all") _pythia.readString("HardQCD:all = on");
    else if (process == "gg") {
      _pythia.readString("HardQCD:gg2gg = on");
      _pythia.readString("HardQCD:qqbar2gg = on");
    }

    _pythia.readString("Tune:pp = "+to_string(tree.get<int>("pythia.tune", 14))); // default tune: Monash
    _pythia.readString("Beams:idA = "+to_string(tree.get<int>("pythia.beamidA", 2212))); // default beam particles: protons
    _pythia.readString("Beams:idB = "+to_string(tree.get<int>("pythia.beamidB", 2212)));
    _pythia.readString("Beams:eCM = "+to_string(tree.get<double>("pythia.eCM", 7000.))); // default center-of-mass energy: 7 TeV
    _pythia.readString("PhaseSpace:pTHatMin = "+to_string(tree.get<double>("pythia.pTHatMin", 100.))); // default minimum pt of generation: 100 GeV
    
    // by default, ISR, FSR, and MPI are all on
    if (!tree.get<bool>("pythia.ISR", true)) _pythia.readString("PartonLevel:ISR = off");
    else _pythia.readString("PartonLevel:ISR = on");
    if (!tree.get<bool>("pythia.FSR", true)) _pythia.readString("PartonLevel:FSR = off");
    else _pythia.readString("PartonLevel:FSR = on");    
    if (!tree.get<bool>("pythia.MPI", true)) _pythia.readString("PartonLevel:MPI = off");
    else _pythia.readString("PartonLevel:MPI = on");

    _is_parton_level = ! tree.get<bool>("pythia.hadronization", true);

    set<string> allowed_hadron_modes = {"all", "measurable", "D/B"};
    _hadron_mode = tree.get<string>("selections.includeHadrons", "D/B");
    if (allowed_hadron_modes.find(_hadron_mode) == allowed_hadron_modes.end()) {
        throw invalid_argument("Invalid mode specified for hadrons: " + _hadron_mode);
    }

    //_include_all_hadrons = tree.get<bool>("pythia.includeAllHadrons", false);

    if (_is_parton_level) {
      _pythia.readString("HadronLevel:all = off");
    }
    else {
      _pythia.readString("HadronLevel:all = on"); // turn on hadronization
      // read through other flags that are specific to hadron-level events, about decays of various hadrons
      // if (_hadron_mode == "D/B") {
      //   if (tree.get<bool>("pythia.switchoffD0decay", true)) _pythia.readString(to_string(constants::D0)+":mayDecay = off");
      //   // by default, B hadrons don't decay:
      //   if (tree.get<bool>("pythia.switchoffbhadrondecays", true)) {
      //     _pythia.readString("511:mayDecay = off");
      //     _pythia.readString("521:mayDecay = off");
      //     _pythia.readString("523:mayDecay = off");
      //     _pythia.readString("513:mayDecay = off");
      //     _pythia.readString("531:mayDecay = off");
      //     _pythia.readString("533:mayDecay = off");
      //     _pythia.readString("5232:mayDecay = off");
      //     _pythia.readString("5122:mayDecay = off");
      //     _pythia.readString("555:mayDecay = off");
      //   }
      // }
      // else if (_hadron_mode == "all") {
      
      for (auto h: chadrons) _pythia.readString(to_string(h)+":mayDecay = off");
      for (auto h: bhadrons) _pythia.readString(to_string(h)+":mayDecay = off");
      
      // }
      if (tree.get<bool>("pythia.nuclearpdfs"), false) {
        _pythia.readString("PDF:useHardNPDFA = on");
        _pythia.readString("PDF:useHardNPDFB = on");
      }
    }

    // by default, pythia generation takes a new seed for each run from the clock time
    if (tree.get<bool>("pythia.random", true)) {
      _pythia.readString("Random:setSeed = on");
      _pythia.readString("Random:seed = 0"); 
    }
    else {
      cout << "is not random?" << endl;
      _pythia.readString("Random:setSeed = on");
      _pythia.readString("Random:seed = 42"); 
    }

    //* parameters for the event selection *//
    // track selections
    _track_cuts.trackPtMin = tree.get<double>("selections.trackPtMin", 0.0);
    _track_cuts.trackEtaCut = tree.get<double>("selections.trackEtaCut", 10.0);
    _track_cuts.HFPtMin = tree.get<double>("selections.HFPtMin", 0.0);
    // jet selections
    _track_cuts.jetR = tree.get<double>("selections.jetR", 0.4); // default jet radius is 0.4
    _track_cuts.JetPtMin = tree.get<double>("selections.jetPtMin", 100.0);
    _track_cuts.JetPtMax = tree.get<double>("selections.jetPtMax", 1.e10);
    _track_cuts.JetEtaMin = tree.get<double>("selections.jetEtaMin", -2.0);
    _track_cuts.JetEtaMax = tree.get<double>("selections.jetEtaMax", 2.0);
    // soft drop parameters
    _track_cuts.zcut = tree.get<double>("selections.zcut", 0.0); // by default, softdrop parameters are zero, so there is no softdrop
    _track_cuts.beta = tree.get<double>("selections.beta", 0.0);

    // jet algorithms
    set<string> allowed_jet_algorithms = {"kt", "antikt", "CA"};
    string jet_algorithm = tree.get<string>("selections.jetAlgorithm", "antikt"); // default jet algorithm: antikt
    if (allowed_jet_algorithms.find(jet_algorithm) == allowed_jet_algorithms.end()) {
        throw invalid_argument("Invalid jet algorithm specified: " + jet_algorithm);
    }
    if (jet_algorithm=="antikt") _jet_algo = fastjet::antikt_algorithm;
    else if (jet_algorithm=="kt") _jet_algo = fastjet::kt_algorithm;
    else if (jet_algorithm=="CA") _jet_algo = fastjet::cambridge_algorithm;

    string jet_recl_algorithm = tree.get<string>("selections.jetReclAlgorithm", "CA"); // default jet reclustering algorithm: CA
    if (allowed_jet_algorithms.find(jet_recl_algorithm) == allowed_jet_algorithms.end()) {
        throw invalid_argument("Invalid jet reclustering algorithm specified: " + jet_recl_algorithm);
    }
    if (jet_recl_algorithm=="antikt") _jet_recl_algo = fastjet::antikt_algorithm;
    else if (jet_recl_algorithm=="kt") _jet_recl_algo = fastjet::kt_algorithm;
    else if (jet_recl_algorithm=="CA") _jet_recl_algo = fastjet::cambridge_algorithm;

    set<string> allowed_FC_modes = {"bare", "standard", "HFradius"};
    _FC_mode = tree.get<string>("selections.FCmode", "bare");
    if (allowed_FC_modes.find(_FC_mode) == allowed_FC_modes.end()) {
        throw invalid_argument("Invalid mode specified for Flavor Cone: " + _FC_mode);
    }
    _match_splitting = tree.get<bool>("selections.matchSplit", false);
    _apply_sd_to_all = tree.get<bool>("selections.applySDtoall", false);
    _recursive_daughters = tree.get<bool>("selections.recursiveDaughters", false);

    //* parameters for the medium modification *//
    bool do_medium_mod = tree.get<bool>("medium.doMediumModification", false);
    if (do_medium_mod) { // read parameters for the medium modification
      _medium_params.qL = tree.get<double>("medium.qhatL", 4.0);
      _medium_params.L = tree.get<double>("medium.L", 4.0) / constants::invGeVtofm; // converts value provided in fm into GeV^-1
      _medium_params.resolutionMode = tree.get<string>("medium.resolution","thetac");
      set<string> allowed_resolution_modes = {"jet", "zero", "thetac"};
      if (allowed_resolution_modes.find(_medium_params.resolutionMode) == allowed_resolution_modes.end()) {
        throw invalid_argument("Invalid mode specified for resolution: " + _medium_params.resolutionMode);
      }

      _do_energy_loss = tree.get<bool>("medium.doEnergyLoss", false);
      if (_do_energy_loss) { // if you are going to do energy loss, read the extra necessary parameters

        if (!_is_parton_level) throw invalid_argument("Invalid parameter selection: energy loss is currently only possible for parton-level events!");

        _medium_params.omega_c = tree.get<double>("medium.omegac", 60.);
        _medium_params.n = tree.get<double>("medium.n", 6.);
        _medium_params.T = tree.get<double>("medium.T", 0.3);
        _medium_params.alpha_med = tree.get<double>("medium.alphaMed", 0.1);
        _medium_params.jetR = _track_cuts.jetR; // jetR is also needed for the medium analysis
      }
    }

    //* parameters for the event selection *//
    set<string> allowed_event_types = {"inc", "inclusive", "cc", "bb", "qq"};
    string event_type = tree.get<string>("selections.eventSelection");
    if (allowed_event_types.find(event_type) == allowed_event_types.end()) {
        throw invalid_argument("Invalid event type specified: " + event_type);
    }

    if (event_type == "inclusive" || event_type=="inc") {
      _is_inclusive = true;
      _file_label = "inc";
    }
    else if (event_type == "qq") {

      if (!_is_parton_level) throw invalid_argument("Invalid parameter selection: selecting qq events with hadronization on are incompatible selections!");

      _is_inclusive = false;
      _parton_ids = {constants::DOWN,constants::ANTIDOWN};
      _particle_ids = {constants::DOWN,constants::ANTIDOWN};
      _medium_params.mc2 = 0.;
      _file_label = "qq";
    }
    else if (event_type == "cc") {
     _is_inclusive = false;
     _parton_ids = {constants::CHARM,constants::ANTICHARM};
     if (_is_parton_level) _particle_ids = {constants::CHARM,constants::ANTICHARM};
     else {
      if (_hadron_mode == "D/B") _particle_ids = {constants::D0,constants::D0BAR};
      else if (_hadron_mode == "measurable") _particle_ids = {constants::D0,constants::D0BAR};
      else if (_hadron_mode == "all") {
        for (int hadronid: chadrons) {
          _particle_ids.push_back(hadronid);
          _particle_ids.push_back(-hadronid);
        }
      }
     }
     _file_label = "cc";
     _medium_params.mc2 = pow(constants::CHARM_MASS, 2.0);      
    }
    else if (event_type == "bb") {
      _is_inclusive = false;
      _parton_ids = {constants::BOTTOM,constants::ANTIBOTTOM};
      if (_is_parton_level) _particle_ids = {constants::BOTTOM,constants::ANTIBOTTOM};
      else {
        if (_hadron_mode == "D/B") _particle_ids = {constants::B0,constants::B0BAR};
        else if (_hadron_mode == "measurable") _particle_ids = {constants::B0,constants::B0BAR};
        else if (_hadron_mode == "all") {
          for (int hadronid: bhadrons) {
            _particle_ids.push_back(hadronid);
            _particle_ids.push_back(-hadronid);
          }
        }
      }
      _file_label = "bb";
      _medium_params.mc2 = pow(constants::BOTTOM_MASS, 2.0); 
    }


  _error_log.open("logfile_"+_file_label+to_string(label));
  assert(_error_log.is_open() && "cannot open the specified error log file.");

  _pythia.init();
}


/**
 * @brief function to declare histograms for the analysis
 * 
 * @param histogram_names: vector of string handles that are used to specify each histogram uniquely
 */
void globalAnalysis::declare_histograms(vector<string> histogram_names) {

  for (auto name: histogram_names) {
    _histograms.insert(pair<string, Histogram>(name, Histogram(_track_cuts.eecMin, _track_cuts.eecMax, _track_cuts.eec_bin_size, name+"_"+_file_label)));
  }
}


/**
 * @brief function to normalize histograms to cross sections (in picobarns)
 */
void globalAnalysis::normalize_histograms() {

  // normalized to cross section 
  double sigma = _pythia.info.sigmaGen()*1.0e9; // total cross section in pb
  double norm = sigma / _n_events / _track_cuts.eec_bin_size;

  for (auto& histogram_pair : _histograms) {
    histogram_pair.second.multiply(norm);
  }
}


/**
 * @brief function to write all histograms associated with the analysis
 */
void globalAnalysis::write_histograms() {

  for (auto& histogram_pair : _histograms) {
    histogram_pair.second.write();
  }
}
