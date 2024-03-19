#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "functions_ccbar.h"
#include <algorithm>
#include <vector>
#include <numeric>
#include "fastjet/ClusterSequence.hh"


using namespace Pythia8;
using namespace std;


MyUserInfo::MyUserInfo(){
  _pdg_id = -99999;
}

MyUserInfo::MyUserInfo(int pdg_id_in, int global_index_in){
  _pdg_id = pdg_id_in;
  _global_index = global_index_in;
}


//_____________________________________________________________________________
//___return true if value is contained in the vector v, and false otherwise____
//_____________________________________________________________________________
template <typename T>
bool is_contained( T value, vector<T>& v) {
	
  return std::find(v.begin(), v.end(), value) != v.end();
}

//_____________________________________________
//___sum all values in a vector of integers____
//_____________________________________________
int total( vector<int>& myvector ) {

  if (myvector.size()>0) return accumulate(myvector.begin(),myvector.end(),0);
  else return 0;
}

//__________________________________________________________________________________
//___ convert a Pythia particle to a pseudojet object (for use e.g. with fastjet)___
//__________________________________________________________________________________
fastjet::PseudoJet as_pseudojet( Pythia8::Particle particle ){
  return fastjet::PseudoJet( particle.px(), particle.py(), particle.pz(), particle.e() );
}

// TODO: maybe there's a way to make a class with this pseudojet type and also the global index type..? and maybe the values as well
splitting as_splitting( Event& event, splitindex myindex ) {
	
  splitting mysplit = empty_splitting();
  
  if (valid_splitting(myindex)) {
    //cout << "has const: " << as_pseudojet(event[ myindex.in ]).perp() << endl;
    //cout << "empty: " << fastjet::PseudoJet().perp() << endl;
    mysplit = { as_pseudojet(event[ myindex.in ]), as_pseudojet(event[ myindex.out1 ]), as_pseudojet(event[ myindex.out2 ])};
  }
  return mysplit;
}

//________________________________________________________________________________
//___returns an empty splitting object. The empty pseudojet has no constituents___
//________________________________________________________________________________
splitting empty_splitting() {
  splitting mysplit {fastjet::PseudoJet(), fastjet::PseudoJet(), fastjet::PseudoJet()};
  return mysplit;
}

//___________________________________________________________________________
//___determine whether a splitindex object refers to a valid splitting_______
//___convention is that the indices 0 or -1 refer to an invalid splitting____
//___________________________________________________________________________
bool valid_splitting( splitindex splitting ) {

  if (splitting.in==0 && splitting.out1==0 && splitting.out2==0 ) return false;
  else if (splitting.in==-1 && splitting.out1==-1 && splitting.out2==-1 ) return false;
  else return true;
}

//______________________________________________________________________________________________
//___overloaded function to determine whether a splitting object refers to a valid splitting____
//___convention is that splittings are empty pseudojets if they are not valid splittings________
//______________________________________________________________________________________________
bool valid_splitting( splitting mysplit ) {

  // if ( !mysplit.out1.has_constituents() || !mysplit.out2.has_constituents() ) return false;
  if ( mysplit.out1.perp()>0.0 && mysplit.out2.perp()>0.0 ) return true;
  else return false;
}


//______________________________________________________________________________________________
//____return true if the jet is within the specified range in pseudorapidity, false otherwise___
//______________________________________________________________________________________________
bool eta_cut(fastjet::PseudoJet jet, double eta_min, double eta_max)
{
   if (jet.eta() > eta_max || jet.eta() < eta_min) return false;
   else return true;
}

//______________________________________________________________________________________________
//____return momentum fraction z between the highest-pt subjet of a jet and the jet itself______
//______________________________________________________________________________________________
double get_z( fastjet::PseudoJet leading_subjet, fastjet::PseudoJet jet) {

  return (jet.px()*leading_subjet.px() + jet.py()*leading_subjet.py() + jet.pz()*leading_subjet.pz() ) / ( jet.px()*jet.px() + jet.py()*jet.py() + jet.pz()*jet.pz() );
}

//______________________________________________________________________________________________
//____return the transverse momentum kt between the highest-pt subjet and the jet itself________
//______________________________________________________________________________________________
double get_kt( fastjet::PseudoJet leading_subjet, fastjet::PseudoJet jet) {
 
  return jet.kt_distance(leading_subjet); // TODO: check that kt_distance gives the correct value
}

//____________________________________________________________________________________________
//____write the attributes associated with the splitting object mysplit to the file outfile___
//____attributes are defined inside splitting_values__________________________________________
//____________________________________________________________________________________________
int write_splitting( splitting mysplit, ofstream& outfile) {

  struct phys values_gcc= {0};
  splitting_values(mysplit, &values_gcc); //virtuality, Eg, _formtime); 

  int noutput = 8;
  double outputnumb[noutput] = {values_gcc.rccbar, values_gcc.kt, values_gcc.z, values_gcc.Eg, values_gcc.gluonpt, values_gcc.virtuality, values_gcc.eta, values_gcc.phi}; 

  outfile << outputnumb[0];
  for (int i = 1; i < noutput; i++) {
    outfile << " " << outputnumb[i];
  } 
  outfile << endl;

  return 0;
}

//____________________________________________________________________________________________
//____calculate the attributes associated with the splitting object and store them in values__
//____________________________________________________________________________________________
int splitting_values( splitting mysplit, struct phys *values) {

  fastjet::PseudoJet jetin, jetout1, jetout2;
  jetin = mysplit.in;
  jetout1 = mysplit.out1;
  jetout2 = mysplit.out2;
    
  if (valid_splitting(mysplit)) {
    values->rccbar =     jetout1.delta_R( jetout2 );
    values->gluonpt =    jetin.pt();
    values->z =          get_z( jetout1, jetin );
    values->virtuality = jetin.m();
    values->kt =         get_kt( jetout1, jetin );
    values->Eg =         jetin.e();
    values->eta =        jetin.eta();
    values->phi =        jetin.phi();
  }

  return 0;
}

//____________________________________________________________________________________________
//____function to determine the instances of particles with pdg ids in pdglist inside the jet_
//____returns [n,m] for jets with n particles with pdg id pdglist[0] and m particles with id__ //____pdglist[1]
//____________________________________________________________________________________________
vector<int> jet_contents_by_pdg( fastjet::PseudoJet& jet, vector<int>& pdglist, vector<fastjet::PseudoJet>& tagged_particles ) {
  
  vector<fastjet::PseudoJet> particles = jet.constituents();
  
  int pdg;
  vector<int> n_pdgs(pdglist.size(), 0);
  
  if (particles.size()==0) return n_pdgs;
  
  for (int ic=0; ic<particles.size(); ic++) {
    
    for (int ip=0; ip<n_pdgs.size(); ip++) {
      
      pdg = pdglist[ip];
      if ( particles[ic].user_info<MyUserInfo>().pdg_id()==pdg ) {
	n_pdgs[ip]+=1;
	if (!tagged_particles[ip].has_constituents()) tagged_particles[ip] = particles[ic];
      }
    }
  }
  return n_pdgs;	
}

//____________________________________________________________________________________________
//____function to determine the instances of particles in the list hadrons inside the jet_____
//____returns [n,m] for jets with n particles with the same global index as hadrons[0] (which_
//____is 0 or 1 since the global index is unique for particles in the event) and m particles__
//____with the same global index as hadrons[1]
//____________________________________________________________________________________________
vector<int> jet_contents_by_index(fastjet::PseudoJet& jet, vector<fastjet::PseudoJet>& hadrons) {

   vector<int> n_pdgs(hadrons.size(), 0);
   vector<fastjet::PseudoJet> constituents = sorted_by_pt(jet.constituents());
   
  if (constituents.size()==0) return n_pdgs;

  for (long unsigned int ic=0; ic<constituents.size(); ic++) {
    
    for (int ip=0; ip<n_pdgs.size(); ip++) {

      int index = hadrons[ip].user_info<MyUserInfo>().global_index();
      if ( constituents[ic].user_info<MyUserInfo>().global_index()==index ) {

        n_pdgs[ip]+=1;
      }
    }
  }
  return n_pdgs;
}

//____________________________________________________________________________________________
//____read all particles from a pythia event and save final state particles (hadrons) in particles_final
//____and final state partons in particles_final_partonlevel.
//____for efficiency, at this stage also save the global indices of particles that have particle ids 
//____of interest for the analysis in tagged_particles and tagged_particles_partonlevel
//____________________________________________________________________________________________
void read_event(Event& event, vector<fastjet::PseudoJet>& particles_final, vector<fastjet::PseudoJet>& particles_partonlevel,vector<fastjet::PseudoJet>& tagged_particles, vector<fastjet::PseudoJet>& tagged_particles_partonlevel, vector<int>& particle_ids, vector<int>& particle_ids_partonlevel, struct trackCuts& track_cuts) {
	
  int particle_id;
  
  // loop over all particles in the event
  for (int i = 0; i < event.size(); ++i) if (event[i].isFinal() || event[i].isFinalPartonLevel())  {		
					     
      fastjet::PseudoJet particle( event[i].px(), event[i].py(),
				   event[i].pz(), event[i].e() );
      
      // perform cuts to remove particles that would have gone out of the detector
      if (event[i].pT() < track_cuts.trackLowPtCut) continue;
      if (abs(event[i].eta()) > track_cuts.trackEtaCut) continue;
      
      // store the particle id as part of the information about the particle
      particle_id = event[i].id();
      MyUserInfo *myUserInfo = new MyUserInfo(particle_id, i);
      particle.set_user_info(myUserInfo);
      
      
      if (event[i].isFinal()) {
	particles_final.push_back(particle);
	
	if ( particle_ids.size()!=0 && is_contained(particle_id, particle_ids) && particle.perp()>track_cuts.HFLowPtCut ) {
	  tagged_particles.push_back(particle);
	}			
      }
      else if (event[i].isFinalPartonLevel()) {
	
	particles_partonlevel.push_back(particle);
	
	if ( particle_ids_partonlevel.size()!=0 && is_contained(particle_id, particle_ids_partonlevel) && particle.perp()>track_cuts.HFLowPtCut ) {
	  tagged_particles_partonlevel.push_back(particle);
	  // cout << "in read: " << i << " (" << particle_id << "), pt = " << particle.perp() << endl;
	}
      }
    } 
}

		
//_________________________________________________________________________________________________
//____iteratively recluster all particles in jet, using the Cambridge/Aachen algorithm.____________
//____perform the reclustering until the point where one of the particles in tagged_list is in one subjet and the other is in the other, indicating the splitting where they were produced together
//____write out the kinematics of that splitting in outfile
//____________________________________________________________________________________________
splitting iterative_reclustering_doubletag(fastjet::PseudoJet& jet, vector<fastjet::PseudoJet>& tagged_list, double jetR) { 
  
  // set the cambridge/aachen algorithm for jet reclustering       
  fastjet::JetAlgorithm jet_reclustering_algorithm(fastjet::cambridge_algorithm);
  // set the recombination scheme for jet reclustering                       
  fastjet::RecombinationScheme recombScheme_reclustering = fastjet::E_scheme; 
  // use a large jet radius for the reclustering to ensure that all particles in the jet are reclustered
  double jetR_reclustering = 5.0 * jetR;     
  fastjet::JetDefinition jet_def_recluster(jet_reclustering_algorithm, jetR_reclustering, recombScheme_reclustering,fastjet::Best);
  vector< fastjet::PseudoJet > emptySplitting = {fastjet::PseudoJet(), fastjet::PseudoJet(), fastjet::PseudoJet()};

                          
  std::vector<fastjet::PseudoJet> jet_constituents = jet.constituents();
  
  fastjet::ClusterSequence reclustering(jet_constituents, jet_def_recluster);
  fastjet::PseudoJet reclustered_jet = (reclustering.inclusive_jets(0.0)).at(0);  // get the reclustered jet                                                        
    
  fastjet::PseudoJet subjet_1;
  fastjet::PseudoJet subjet_2;
  vector<int> contents_1, contents_2;
    
  while (reclustered_jet.has_parents( subjet_1, subjet_2)) {
        
    // get the one-hot encoding for the particles in tagged_list in each subjet; namely, if the particle tagged_list[0] is contained in subjet_1 then contents_1[0] will be 1, otherwise 0
    contents_1 = jet_contents_by_index( subjet_1, tagged_list );
    contents_2 = jet_contents_by_index( subjet_2, tagged_list );
    
    // if each subjet has exactly one of the particles from tagged_list (D and Dbar, or c and cbar), and between them they have both, then break out of the loop and return this splitting
    if ( total(contents_1)==1 && total(contents_2)==1 && ( (contents_1[0]==1&&contents_2[1]==1) || (contents_1[1]==1&&contents_2[0]==1) ) ) 
      break;

    // if both Ds are in subjet 2, then redefine so that they are in subjet 1 (since we will follow that one)
    if ( total(contents_2)==2 && total(contents_1)==0 ) 
      swap( subjet_1, subjet_2 );

    reclustered_jet = subjet_1; // reiterate with the subjet that contains both particles                                                  
  }

  // now that each subjet contains on of the particles in tagged_list, define subjet 1 to be the higher pt one
  if ( subjet_2.perp()>subjet_1.perp() ) swap( subjet_1, subjet_2 );

  return {reclustered_jet, subjet_1, subjet_2};
}


//_________________________________________________________________________________________________
//____function to initialize pythia given the parameters for the run given in the parameter file___
//____also sets values for elements of trackCuts struct based on the parameter file________________
//____turns off decays of a lot of B hadrons as well, to make sure they don't decay into Ds________
//____returns the number of events, for use in the main function___________________________________
//_________________________________________________________________________________________________
int initialize_pythia( string paramfile_name, Pythia& pythia, struct trackCuts& track_cuts) {
	
  fstream paramfile;
  paramfile.open(paramfile_name, ios::in);
  string param_name, param_value;
  int maxnevents=0;
  bool switchoffbhadrondecays;

  // process name
  pythia.readString("HardQCD:all = on");
  
  // initialize pythia according to the provided values in paramfile
  while ( paramfile >> param_name >> param_value )
    {
      if (param_name == "maxneventsperjob:") {
	maxnevents = stoi(param_value);
	pythia.readString("Main:numberOfEvents = "+param_value);
      }
      else if (param_name == "tune:") pythia.readString("Tune:pp = "+param_value);
      else if (param_name == "beamidA:") pythia.readString("Beams:idA = "+param_value);
      else if (param_name == "beamidB:") pythia.readString("Beams:idB = "+param_value);
      else if (param_name == "eCM:") pythia.readString("Beams:eCM = "+param_value);
      else if (param_name == "pTHatMin:") pythia.readString("PhaseSpace:pTHatMin = "+param_value);
      else if (param_name == "jetR:") track_cuts.jetR = stof(param_value);
      else if (param_name == "trackLowPtCut:") track_cuts.trackLowPtCut = stof(param_value);
      else if (param_name == "trackEtaCut:") track_cuts.trackEtaCut = stof(param_value);
      else if (param_name == "HFLowPtCut:") track_cuts.HFLowPtCut = stof(param_value);
      else if (param_name == "pThat:") pythia.readString("PhaseSpace:pTHatMin = "+param_value);
      else if (param_name == "jetminpt:") track_cuts.JetLowPtCut = stof(param_value);
      else if (param_name == "ISR:") pythia.readString("PartonLevel:ISR = "+param_value);
      else if (param_name == "FSR:") pythia.readString("PartonLevel:FSR = "+param_value);	
      else if (param_name == "switchoffbhadrondecays:") {
	if (param_value=="true") switchoffbhadrondecays = true;
	else if (param_value=="false") switchoffbhadrondecays = false;
	else cout << "invalid parameter: switchoffbhadrondecays should be true or false" << endl;		
      }
 }
  
  // turn off decays of some hadrons, partiuclarly the D0 and various B hadrons that could decay to D's
  pythia.readString(to_string(constants::D0)+":mayDecay = off");
  if (switchoffbhadrondecays) {
    pythia.readString("511:mayDecay = off");
    pythia.readString("521:mayDecay = off");
    pythia.readString("523:mayDecay = off");
    pythia.readString("513:mayDecay = off");
    pythia.readString("531:mayDecay = off");
    pythia.readString("533:mayDecay = off");
    pythia.readString("5232:mayDecay = off");
    pythia.readString("5122:mayDecay = off");
    pythia.readString("555:mayDecay = off");
  }
	
  pythia.init();
  return maxnevents;
}


//_________________________________________________________________________________________________
//____function to go through all jets, and do iterative declustering on any jet containing one of__
//____each of the particles contained in tagged_particles (in standard usage, D/Dbar pair or_______
//____c/cbar pair). Then writes the properties of thesplitting between the D/Dbar (or c/cbar) pair_ //____found by iterative declustering into outfile_________________________________________________
//_________________________________________________________________________________________________
void process_jets(vector<fastjet::PseudoJet>& final_particles, vector<fastjet::PseudoJet>& tagged_particles, vector<int>& tagged_particle_ids, ofstream& outfile, struct trackCuts& track_cuts, int iEvent ) {
	
  vector<fastjet::PseudoJet> jets;
  fastjet::PseudoJet jet;
	
  // to avoid doing extra work, reject any cases where there are not at least two tagged particles in the event, since we will eventually look for jets containing a D/Dbar (c/cbar) pair
  if ( tagged_particles.size()<2 ) return;

  // cluster jets	using the anti-kt algorithm with minimum pt given by jetminpt
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, track_cuts.jetR );	
  fastjet::ClusterSequence cluster_seq(final_particles, jet_def);
  jets = sorted_by_pt( cluster_seq.inclusive_jets(track_cuts.JetLowPtCut) );
	
  // consider at most the two highest-pt jets in the event
  int njets=2;
  if (jets.size()<2) njets=jets.size();
  vector<int> one_each{1,1};
	
  // loop over those jets
  for (int ijet = 0; ijet < njets; ijet++) {

    jet = jets.at(ijet);
 
    // reject jets that do not fall in the pseudorapidity acceptance cuts
    if (!eta_cut(jet, track_cuts.JetEtaMin, track_cuts.JetEtaMax) ) continue;

    // reject jets that have no constituents
    if ( !jet.has_constituents() ) continue;

    // reject jets that do not have at least one particle associated with each of the particle ids provided
    if ( !(jet_contents_by_pdg( jet, tagged_particle_ids, tagged_particles )==one_each) ) continue;

    splitting tagged_split = iterative_reclustering_doubletag( jet, tagged_particles, track_cuts.jetR );

    outfile << iEvent << " ";
    write_splitting( tagged_split, outfile);
  
  }
}


//_________________________________________________________________________________________________
//____function to iterate through the tree structure provided by the Monte Carlo generator, to_____
//____obtain approximate "truth level" information to compare to the reclustering on final-state___
//____particles only. Writes out the information for the splitting that produced tagged_partons,___
//____according to the Monte Carlo generator itself________________________________________________
//_________________________________________________________________________________________________
void process_MC(Event& event, vector<fastjet::PseudoJet> tagged_partons, ofstream& outfile, int iEvent) {

      splitting mysplit = find_common_splitting(event, tagged_partons[0].user_info<MyUserInfo>().global_index(), tagged_partons[1].user_info<MyUserInfo>().global_index());

      outfile << iEvent << " ";
      write_splitting( mysplit, outfile );

}

//_________________________________________________________________________________________________
//____return the index of the sister of the particle indexed by particle_index in the event. The___
//____sister is the particle that has the same mother as the particle._____________________________
//_________________________________________________________________________________________________
int sister(Event& event, int particle_index) {
  
  int mother_index, daugh1_index, daugh2_index;
  
  mother_index = event[particle_index].mother1();
  daugh1_index = event[ mother_index ].daughter1();
  daugh2_index = event[ mother_index ].daughter2();

  // whichever one is not the current particle, is the sister
  if (daugh1_index==particle_index) return daugh2_index;
  else if (daugh2_index==particle_index) return daugh1_index;
  else return -1;
}

//_________________________________________________________________________________________________
//____return the splitting which is the most recent splitting that produced particles 1 and 2,_____
//____namely, the most recent shared ancestor between particles 1 and 2.___________________________
//_________________________________________________________________________________________________
splitting find_common_splitting(Event& event, int particle_1, int particle_2) {

  int current_particle, current_mother;
  vector<int> incoming_status {constants::INCOMING_HARD, constants::INCOMING_SUB, constants::INCOMING_ISR};
  
  current_particle = particle_1;
  current_mother = event[ current_particle ].mother1();

  // while the current particle is not an ancestor of the current particle, and the current particle does not have a status code associated with being an incoming particle
  while ( !event[ particle_2 ].isAncestor( current_particle ) && !is_contained( event[ current_particle ].status(), incoming_status ) ) {
  
    current_particle = event[ current_particle ].mother1();
  }

  splitindex split = { current_particle, event[current_particle].daughter1(), event[current_particle].daughter2() };

  return as_splitting(event, split);
}

//_________________________________________________________________________________________________
//____return the iterator in a vector of PseudoJet objects corresponding to the first instance_____
//____of the pdg id of the PseudoJet having the value given as an argument.________________________
//_________________________________________________________________________________________________
vector<fastjet::PseudoJet>::iterator find_first_instance(vector<fastjet::PseudoJet>& vec, const int& value) {
    for (auto it = vec.begin(); it != vec.end(); ++it) {
      if ((*it).user_info<MyUserInfo>().pdg_id() == value) {
            return it;
        }
    }
    return vec.end(); // Return end iterator if no such element is found
}

//_________________________________________________________________________________________________
//____return the highest-pt pair from tagged particles that has the highest pt of each of the______
//____particles indexed by tagged_pids. For example, this gives the highest pt c/cbar (D/Dbar) pair
//_________________________________________________________________________________________________
vector<fastjet::PseudoJet> highest_pt_pair( vector<fastjet::PseudoJet> tagged_particles, vector<int> tagged_pids ) {

  vector<fastjet::PseudoJet> sorted_tagged = sorted_by_pt(tagged_particles);

  vector<fastjet::PseudoJet> objs;
  vector<fastjet::PseudoJet>::iterator it;
  for (int i=0; i<tagged_pids.size(); i++) {

    it = find_first_instance( sorted_tagged, tagged_pids[i] );
    if (it!=sorted_tagged.end()) objs.push_back(*it);
  }
  
  return objs;
}
