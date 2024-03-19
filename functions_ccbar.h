#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"

using namespace Pythia8;
using namespace std;


// structures
struct phys {
  double Eg;
  double virtuality;
  double z;
  double rccbar;
  double gluonpt;
  double kt;
  double eta;
  double phi;
};

struct splitindex { // object containing global (event-wide) indices of the splitting candidates
  int in, out1, out2;
};

struct splitting { // object containing PseudoJet objects of the splitting candidates
  fastjet::PseudoJet in, out1, out2;
};

// object containing track, jet, and heavy flavor cuts and acceptances for all events
struct trackCuts {
  double trackEtaCut;
  double trackLowPtCut;
  double jetR;
  double HFLowPtCut;
  double JetLowPtCut;
  double JetEtaMin;
  double JetEtaMax;
};

//------------------------------------------------------------------------
// the user class
//
// To associate extra information to a PseudoJet, one first has to
// create a class, derived from UserInfoBase, that contains
// that information.
//
// To associate a user-defined extra information, PseudoJet::set_user_info()
//
class MyUserInfo : public fastjet::PseudoJet::UserInfoBase
{
 public:
  //  pdg_id: the PDG id of the particle
  //  global_index: the index of the particle in the event record
  MyUserInfo();
  MyUserInfo(int pdg_idin, int global_indexin);

  /// access to the PDG id
  int pdg_id() const { return _pdg_id; }

  // access to the event-level index
  int global_index() const { return _global_index; }

  int _pdg_id;       // the associated pdg id
  int _global_index; // the associated event index
};

namespace constants
{
  // identifying particle ids used by pythia
  const int CHARM = 4;
  const int ANTICHARM = -4;
  const int GLUON = 21;
  const int D0 = 421;
  const int D0BAR = -421;
  // status code for particles in pythia
  const int INCOMING_HARD = -21; // incoming from hard process
  const int INCOMING_SUB = -31; // incoming from subprocess
  const int INCOMING_ISR = -41; // incoming from initial state radiation
  const int SHOWER = -51; // produced by the shower

}

//_____________________________________________________________________________
//___return true if value is contained in the vector v, and false otherwise____
//_____________________________________________________________________________
template <typename T>
bool is_contained( T value, vector<T>& v);

//_____________________________________________
//___sum all values in a vector of integers____
//_____________________________________________
int total( vector<int>& myvector );

//__________________________________________________________________________________
//___ convert a Pythia particle to a pseudojet object (for use e.g. with fastjet____
//__________________________________________________________________________________
fastjet::PseudoJet as_pseudojet( Pythia8::Particle particle );

splitting as_splitting( Event& event, splitindex myindex );

//________________________________________________________________________________
//___returns an empty splitting object. The empty pseudojet has no constituents___
//________________________________________________________________________________
splitting empty_splitting();

//___________________________________________________________________________
//___determine whether a splitindex object refers to a valid splitting_______
//___convention is that the indices 0 or -1 refer to an invalid splitting____
//___________________________________________________________________________
bool valid_splitting( splitindex splitting );

//______________________________________________________________________________________________
//___overloaded function to determine whether a splitting object refers to a valid splitting____
//___convention is that splittings are empty pseudojets if they are not valid splittings________
//______________________________________________________________________________________________
bool valid_splitting( splitting mysplit );

//______________________________________________________________________________________________
//____return true if the jet is within the specified range in pseudorapidity, false otherwise___
//______________________________________________________________________________________________
bool eta_cut(fastjet::PseudoJet jet, double eta_min, double eta_max);

//______________________________________________________________________________________________
//____return momentum fraction z between the highest-pt subjet of a jet and the jet itself______
//______________________________________________________________________________________________
double get_z( fastjet::PseudoJet leading_subjet, fastjet::PseudoJet jet);

//______________________________________________________________________________________________
//____return the transverse momentum kt between the highest-pt subjet and the jet itself________
//______________________________________________________________________________________________
double get_kt( fastjet::PseudoJet leading_subjet, fastjet::PseudoJet jet);

//____________________________________________________________________________________________
//____write the attributes associated with the splitting object mysplit to the file outfile___
//____attributes are defined inside splitting_values__________________________________________
//____________________________________________________________________________________________
int write_splitting( splitting mysplit, ofstream& outfile);

//____________________________________________________________________________________________
//____calculate the attributes associated with the splitting object and store them in values__
//____________________________________________________________________________________________
int splitting_values( splitting mysplit, struct phys *values);

//____________________________________________________________________________________________
//____function to determine the instances of particles with pdg ids in pdglist inside the jet_
//____returns [n,m] for jets with n particles with pdg id pdglist[0] and m particles with id__ //____pdglist[1]
//____________________________________________________________________________________________
vector<int> jet_contents_by_pdg( fastjet::PseudoJet& jet, vector<int>& pdglist, vector<fastjet::PseudoJet>& tagged_particles );

//____________________________________________________________________________________________
//____function to determine the instances of particles in the list hadrons inside the jet_____
//____returns [n,m] for jets with n particles with the same global index as hadrons[0] (which_
//____is 0 or 1 since the global index is unique for particles in the event) and m particles__
//____with the same global index as hadrons[1]
//____________________________________________________________________________________________
vector<int> jet_contents_by_index( fastjet::PseudoJet& jet, vector<fastjet::PseudoJet>& hadrons );

//____________________________________________________________________________________________
//____read all particles from a pythia event and save final state particles (hadrons) in particles_final
//____and final state partons in particles_final_partonlevel.
//____for efficiency, at this stage also save the global indices of particles that have particle ids 
//____of interest for the analysis in tagged_particles and tagged_particles_partonlevel
//____________________________________________________________________________________________
void read_event(Event& event, vector<fastjet::PseudoJet>& particles_final, vector<fastjet::PseudoJet>& particles_partonlevel,vector<fastjet::PseudoJet>& tagged_particles, vector<fastjet::PseudoJet>& tagged_particles_partonlevel, vector<int>& particle_ids, vector<int>& particle_ids_partonlevel, struct trackCuts& track_cuts);

//_________________________________________________________________________________________________
//____iteratively recluster all particles in jet, using the Cambridge/Aachen algorithm.____________
//____perform the reclustering until the point where one of the particles in tagged_list is in one subjet and the other is in the other, indicating the splitting where they were produced together
//____write out the kinematics of that splitting in outfile
//____________________________________________________________________________________________
splitting iterative_reclustering_doubletag(fastjet::PseudoJet& jet, vector<fastjet::PseudoJet>& tagged_list, double jetR);

//_________________________________________________________________________________________________
//____function to initialize pythia given the parameters for the run given in the parameter file___
//____also sets values for elements of trackCuts struct based on the parameter file________________
//____turns off decays of a lot of B hadrons as well, to make sure they don't decay into Ds________
//_________________________________________________________________________________________________
int initialize_pythia( string paramfile_name, Pythia& pythia, struct trackCuts& track_cuts);

//_________________________________________________________________________________________________
//____function to go through all jets, and do iterative declustering on any jet containing one of__
//____each of the particles contained in tagged_particles (in standard usage, D/Dbar pair or_______
//____c/cbar pair). Then writes the properties of thesplitting between the D/Dbar (or c/cbar) pair_ //____found by iterative declustering into outfile_________________________________________________
//_________________________________________________________________________________________________
void process_jets(vector<fastjet::PseudoJet>& final_particles, vector<fastjet::PseudoJet>& tagged_particles, vector<int>& tagged_particle_ids, ofstream& outfile, struct trackCuts& track_cuts, int iEvent );

//_________________________________________________________________________________________________
//____function to iterate through the tree structure provided by the Monte Carlo generator, to_____
//____obtain approximate "truth level" information to compare to the reclustering on final-state___
//____particles only. Writes out the information for the splitting that produced tagged_partons,___
//____according to the Monte Carlo generator itself________________________________________________
//_________________________________________________________________________________________________
void process_MC(Event& event, vector<fastjet::PseudoJet> tagged_partons, ofstream& outfile, int iEvent);

//_________________________________________________________________________________________________
//____return the index of the sister of the particle indexed by particle_index in the event. The___
//____sister is the particle that has the same mother as the particle._____________________________
//_________________________________________________________________________________________________
int sister(Event& event, int particle_index);

//_________________________________________________________________________________________________
//____return the splitting which is the most recent splitting that produced particles 1 and 2,_____
//____namely, the most recent shared ancestor between particles 1 and 2.___________________________
//_________________________________________________________________________________________________
splitting find_common_splitting(Event& event, int particle1_index, int particle2_index);

//_________________________________________________________________________________________________
//____return the iterator in a vector of PseudoJet objects corresponding to the first instance_____
//____of the pdg id of the PseudoJet having the value given as an argument.________________________
//_________________________________________________________________________________________________
vector<fastjet::PseudoJet>::iterator find_first_instance(vector<fastjet::PseudoJet>& vec, const int& value);

//_________________________________________________________________________________________________
//____return the highest-pt pair from tagged particles that has the highest pt of each of the______
//____particles indexed by tagged_pids. For example, this gives the highest pt c/cbar (D/Dbar) pair
//_________________________________________________________________________________________________
vector<fastjet::PseudoJet> highest_pt_pair( vector<fastjet::PseudoJet> tagged_particles, vector<int> tagged_pids );
