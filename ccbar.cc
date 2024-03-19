#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "functions_ccbar.h"
#include "fastjet/ClusterSequence.hh"

using namespace Pythia8;
using namespace std;



int main(int argc, char* argv[])
{
	
  string mycase = string("splitting");

  string paramfile_name = string("params.dat");
  Pythia pythia;
  struct trackCuts track_cuts = {0};
  
  int maxnevents = initialize_pythia(paramfile_name, pythia, track_cuts);

  std::cout << "Using R = " << track_cuts.jetR << std::endl;
  ofstream hadronfile(mycase + "_hadron.dat");
  if (!hadronfile.is_open())
    cout << "Unable to open file";

  ofstream partonfile(mycase + "_parton.dat");
  if (!partonfile.is_open())
    cout << "Unable to open file";

    ofstream MCfile(mycase + "_MC.dat");
  if (!MCfile.is_open())
    cout << "Unable to open file";

  // define cuts for tracks, jets, and heavy flavor particles
  track_cuts.JetEtaMin = -track_cuts.trackEtaCut + track_cuts.jetR; // jet eta range
  track_cuts.JetEtaMax = -track_cuts.JetEtaMin;
  track_cuts.HFLowPtCut = 5; // only tag heavy-flavor particles with pT > 5 GeV

  vector<fastjet::PseudoJet> final_hadrons, final_partons; // only final state particles (at hadron and parton-level, respectively)
  vector<fastjet::PseudoJet> tagged_hadrons, tagged_partons;
  vector<int> hadron_ids{constants::D0, constants::D0BAR};
  vector<int> parton_ids{constants::CHARM, constants::ANTICHARM};

  // Begin event loop. Generate event. Skip if error. List first one. 
  for (int iEvent = 0; iEvent < maxnevents; ++iEvent) {
    	  
    // ignore events where pythia aborts
    if (!pythia.next())
      continue;
    Event& event = pythia.event;

    // initialize per event
    final_hadrons.clear();
    final_partons.clear();	 
    tagged_hadrons.clear();
    tagged_partons.clear();
    
	// read in particles from the event. Save final-state hadrons, final-state partons, and "tagged_partons" and "tagged hadrons" corresponding to particles with the particle ids contained in hadron_ids and parton_ids
    read_event( event, final_hadrons, final_partons, tagged_hadrons, tagged_partons, hadron_ids, parton_ids, track_cuts );
    
    if (tagged_hadrons.size()>=2) {

      tagged_hadrons = highest_pt_pair(tagged_hadrons, hadron_ids);
      if (tagged_hadrons.size()!=2) continue;
      
      // process jets associated with final hadrons, and write splittings from the iterative reclustering
      process_jets( final_hadrons, tagged_hadrons, hadron_ids, hadronfile, track_cuts, iEvent);
    }

    if (tagged_partons.size()>=2) {

      tagged_partons = highest_pt_pair(tagged_partons, parton_ids);
      if (tagged_partons.size()!=2) continue;
      
      // same as above except with final partons
      process_jets( final_partons, tagged_partons, parton_ids, partonfile, track_cuts, iEvent);

      process_MC(event, tagged_partons, MCfile, iEvent);
      
    }
  } // End of event loop.

  return 0;
}
