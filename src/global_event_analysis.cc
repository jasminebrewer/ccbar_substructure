#include <fstream> 

#include "global_event_analysis.hh"
#include "constants.hh"


/**
 * @brief function to set parameters of the analysis and initialize pythia from values specified in the parameter file
 */
void globalAnalysis::initialize_pythia() {
	
  string param_name, param_value;
  bool switchoffbhadrondecays=true;
  bool D0decay=false;

  std::string preamble, event_type;

 // initialize pythia according to the provided values in paramfile
  while ( _parameter_file >> param_name >> param_value )
    {
      // save parameters to prepend to log file
      preamble.append( "# " + param_name + param_value + '\n' );
      if (param_name == "maxneventsperjob:") {
	      _n_events = stoi(param_value);
	      _pythia.readString("Main:numberOfEvents = "+param_value);
      }
      else if (param_name == "process:") {
        if (param_value == "all") _pythia.readString("HardQCD:all = on");
        else if (param_value == "gg") {
          _pythia.readString("HardQCD:gg2gg = on");
          _pythia.readString("HardQCD:qqbar2gg = on");
        }
        else cout << "invalid parameter: process can be 'all' (HardQCDall=on) or 'gg' (gluon-only initial states)." << endl;
      }
      else if (param_name == "eventSelection:") event_type = param_value;
      // parameters for physical features in pythia
      else if (param_name == "tune:") _pythia.readString("Tune:pp = "+param_value);
      else if (param_name == "beamidA:") _pythia.readString("Beams:idA = "+param_value);
      else if (param_name == "beamidB:") _pythia.readString("Beams:idB = "+param_value);
      else if (param_name == "eCM:") _pythia.readString("Beams:eCM = "+param_value);
      else if (param_name == "pTHatMin:") _pythia.readString("PhaseSpace:pTHatMin = "+param_value);
      else if (param_name == "ISR:") _pythia.readString("PartonLevel:ISR = "+param_value);
      else if (param_name == "MPI:") _pythia.readString("PartonLevel:MPI = "+param_value);
      else if (param_name == "FSR:") _pythia.readString("PartonLevel:FSR = "+param_value);	
      // parameters to specify analysis cuts
      else if (param_name == "jetR:") _track_cuts.jetR = stof(param_value);
      else if (param_name == "trackPtMin:") _track_cuts.trackPtMin = stof(param_value);
      else if (param_name == "trackEtaCut:") {
        _track_cuts.trackEtaCut = stof(param_value);
      }
      else if (param_name == "HFPtMin:") _track_cuts.HFPtMin = stof(param_value);
      else if (param_name == "jetPtMin:") _track_cuts.JetPtMin = stof(param_value);
      else if (param_name == "jetPtMax:") _track_cuts.JetPtMax = stof(param_value);
      else if (param_name == "jetEtaMin:") _track_cuts.JetEtaMin = stof(param_value);
      else if (param_name == "jetEtaMax:") _track_cuts.JetEtaMax = stof(param_value);
      // hadronization and hadron decays
      else if (param_name == "hadronization:") {
        if (param_value=="off") {
          _pythia.readString("HadronLevel:all = off");
          _is_parton_level = true; }
        else if (param_value=="on") {
          _pythia.readString("HadronLevel:all = on");
          _is_parton_level = false; }
      }
      else if (param_name == "switchoffbhadrondecays:") {
	if (param_value=="true") switchoffbhadrondecays = true;
	else if (param_value=="false") switchoffbhadrondecays = false;
	else cout << "invalid parameter: switchoffbhadrondecays should be true or false" << endl;		
      }
      else if (param_name == "D0decays:") {
	if (param_value=="true") D0decay = true;
	else if (param_value=="false") D0decay = false;
	else cout << "invalid parameter: D0decays should be true or false" << endl;		
      }
      // medium parameters
      else if (param_name == "qhatL:") _qhatL = stof(param_value);
      else if (param_name == "L:") _L = stof(param_value) / constants::invGeVtofm;
      else if (param_name == "jetAlgorithm:") {
        if (param_value=="antikt") _jet_algo = fastjet::antikt_algorithm;
        else if (param_value=="kt") _jet_algo = fastjet::kt_algorithm;
        else if (param_value=="CA") _jet_algo = fastjet::cambridge_algorithm;
      }
      else if (param_name == "jetReclAlgorithm:") {
        if (param_value=="antikt") _jet_recl_algo = fastjet::antikt_algorithm;
        else if (param_value=="kt") _jet_recl_algo = fastjet::kt_algorithm;
        else if (param_value=="CA") _jet_recl_algo = fastjet::cambridge_algorithm;       
      }
      else {
        cout << "unknown command line parameter " << param_name << endl;
      }
    }

    if (event_type == "inclusive") {
      _is_inclusive = true;
      _file_label = "inc";
    }
    else if (event_type == "qq") {
      _is_inclusive = false;
      if (_is_parton_level) _particle_ids = {constants::DOWN,constants::ANTIDOWN};
      else cout << "WARNING: hadronization: on and eventSelection: qq are incompatible selections!" << endl;
      _mc2 = 0.;
      _file_label = "qq";
    }
    else if (event_type == "cc") {
     _is_inclusive = false;
     _parton_ids = {constants::CHARM,constants::ANTICHARM};
     if (_is_parton_level) _particle_ids = {constants::CHARM,constants::ANTICHARM};
     else _particle_ids = {constants::D0,constants::D0BAR};
     _file_label = "cc";
     _mc2 = pow(constants::CHARM_MASS, 2.0);        
    }
    else if (event_type == "bb") {
      _is_inclusive = false;
      _parton_ids = {constants::BOTTOM,constants::ANTIBOTTOM};
      if (_is_parton_level) _particle_ids = {constants::BOTTOM,constants::ANTIBOTTOM};
      else _particle_ids = {constants::B0,constants::B0BAR};
      _file_label = "bb";
      _mc2 = pow(constants::BOTTOM_MASS, 2.0); 
    }
    else { cout << "parameter eventSelection has the allowed values 'inclusive' (all jets), 'qq' (down-quark tagged jets), 'cc' (charm-tagged jets), or 'bb' (bottom-tagged jets)."; }


  _error_log.open("logfile_"+_file_label);
  assert(_error_log.is_open() && "cannot open the specified error log file.");
  _error_log << preamble;
 

  // Pick new random number seed for each run, based on clock
  _pythia.readString("Random:setSeed = on");
  _pythia.readString("Random:seed = 0"); 

  // turn off decays of some hadrons, partiuclarly the D0 and various B hadrons that could decay to D's
  if (!D0decay) _pythia.readString(to_string(constants::D0)+":mayDecay = off");
  if (switchoffbhadrondecays) {
    _pythia.readString("511:mayDecay = off");
    _pythia.readString("521:mayDecay = off");
    _pythia.readString("523:mayDecay = off");
    _pythia.readString("513:mayDecay = off");
    _pythia.readString("531:mayDecay = off");
    _pythia.readString("533:mayDecay = off");
    _pythia.readString("5232:mayDecay = off");
    _pythia.readString("5122:mayDecay = off");
    _pythia.readString("555:mayDecay = off");
  }
	
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
