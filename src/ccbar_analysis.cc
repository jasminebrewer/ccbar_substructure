#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "ccbar_analysis.hh"
#include "histograms.hh"
#include "fastjet/ClusterSequence.hh"


using namespace Pythia8;
using namespace fastjet;
using namespace std;

/**
 * @brief helper function to determine if a value is contained within a vector
 * @param value
 * @param vector: vector with the same type as values
 * @return true if value is contained in vector, false otherwise
 */
template <typename T>
bool is_contained( T value, vector<T>& v) { return std::find(v.begin(), v.end(), value) != v.end(); }


/**
 * @brief helper function to convert a Pythia particle object into a PseudoJet object
 * 
 * converts a pythia particle into a PseudoJet, including setting user_info of the particle id and the global index of the particle in the pythia event
 * 
 * @param particle: pythia particle to convert to psuedojet
 * @return PseudoJet object with the momentum, particle id, and index of particle
 */
PseudoJet as_pseudojet( Particle& particle ){
  PseudoJet pj = PseudoJet( particle.px(), particle.py(), particle.pz(), particle.e() );
  pj.set_user_info(new ExtraInfo(particle.id(), particle.index()));
  return pj;
}


/**
 * @brief helper function to determine if a particle is among the constituents of a jet
 *  
 * a particle is deemed to be among the jet constituents if it shares its global event index with one of the jet constituents
 * 
 * @param jet PseudoJet object clustering one or more particles
 * @param particle PseudoJet object for a particular particle
 * @return true if particle is within the jet, false otherwise
 */
bool contains(PseudoJet& jet, PseudoJet particle) {

  int global_index = particle.user_info<ExtraInfo>().global_index();
  vector<PseudoJet> constituents = sorted_by_pt(jet.constituents());
  if (constituents.size()==0) return false;
  for (auto c: constituents) {
    if (c.user_info<ExtraInfo>().global_index() == global_index ) return true;
  }
  return false;
}

/**
 * @brief function to trace a non-final state pythia particle to the corresponding final state particle
 *  
 * the particle is assumed to have a unique daughter with the same particle id which allows that particle to be traced to the final state
 * (namely, this function should be used for particles undergoing q->qg splittings, and the function follows the q)
 * 
 * @param particle PseudoJet object for a particular particle to trace to the final state
 * @return PseudoJet object associated with the final state version of that particle
 */
PseudoJet EventCCbar::follow_to_final_state( PseudoJet particle) {

  int current_index = particle.user_info<ExtraInfo>().global_index(); // event index of the particle
  int pid = particle.user_info<ExtraInfo>().pdg_id(); // its particle id
  Particle& current_particle = _py_event[ current_index ]; // pythia particle associated with that particle
  
  // while the current particle is not final
  while (!current_particle.isFinal()) {

    current_particle = _py_event[ current_index ];
 
    // set the new current particle to which daughter has the same particle id as the current particle
    if ( _py_event[ current_particle.daughter1() ].id()==pid) current_index = current_particle.daughter1();
    else if ( _py_event[ current_particle.daughter2() ].id()==pid) current_index = current_particle.daughter2();
  }

  // convert the final particle to a PseudoJet object and return
  return as_pseudojet( _py_event[ current_index ] );
}

/**
 * @brief read in a pythia event
 *  
 * from the pythia event, store final state particles passing the track selection cuts, and store global indices of particles of particular interest to the analysis, as class members
 */
void EventCCbar::read_event() {

  double max_pt_particle=0.;
  double max_pt_antiparticle=0.;
  bool has_particle=false;
  bool has_antiparticle=false;

  // loop over all particles in the event
  for (Particle& p: _py_event) {

    if (!p.isFinal()) continue; // consider only final states

    // perform cuts to remove particles that would have gone out of the detector
    if (p.pT() < _track_cuts.trackPtMin) continue;
    if (abs(p.eta()) > _track_cuts.trackEtaCut) continue;
    
    PseudoJet particle( p.px(), p.py(), p.pz(), p.e() );
    
    // store the particle id and global index as part of about the particle
    particle.set_user_info(new ExtraInfo(p.id(), p.index()));
    
    _final_particles.push_back(particle);

    if (_is_inclusive) continue; // no extra particle info to save, so we're done

    // add a particle to tagged_particles if it has a particle id matching particle_ids, and if its pT is bigger than the low pT cut for heavy flavor 
    if (is_contained(p.id(), _particle_ids) && particle.perp()>_track_cuts.HFPtMin ) {
      _tagged_particles.push_back(particle);

      // keep track of the highest-pT tagged particle and tagged antiparticle in the event
      if (p.id()>0 && particle.perp()>max_pt_particle) {
        max_pt_particle = particle.perp();
        _maxpt_tagged_particle = particle;
        has_particle=true;
      }
      if (p.id()<0 && particle.perp()>max_pt_antiparticle) { 
        max_pt_antiparticle = particle.perp();
        _maxpt_tagged_antiparticle = particle;
        has_antiparticle=true;
      }
    }
  }
  _has_pair = has_particle && has_antiparticle;
  // sort the tagged particles by pT
  _tagged_particles = sorted_by_pt(_tagged_particles); 
}


/**
 * @brief cluster all final-state particles in the event into jets
 *  
 * cluster jets according to the provided jet definition, reject jets that don't satisfy the jet selection criteria, and store up to two highest pt jets as class members
 */void EventCCbar::cluster_jets() {

   // cluster final particles in the event according to the jet definition jet_def and save the associated cluster sequence for later
   vector<PseudoJet> all_jets;
   ClusterSequence cs(_final_particles, _jet_def);
   _cluster_seq = cs;
   // store all jets in the event above the min pT cut and sort them by pT
   all_jets = sorted_by_pt( _cluster_seq.inclusive_jets(_track_cuts.JetPtMin) );

   for (auto jet: all_jets) {

     // reject jets that either don't have constituents or fall outside of the acceptance cuts of the analysis
     bool pass_jet_cuts = (jet.has_constituents() && jet.eta()>_track_cuts.JetEtaMin && jet.eta()<_track_cuts.JetEtaMax && jet.pt()>_track_cuts.JetPtMin && jet.pt()<_track_cuts.JetPtMax );
     if (!pass_jet_cuts) continue;

     // add at most 2 highest pt jets from the event satisfying the cuts to the event's jets
     if ( _jets.size() < 2 ) _jets.push_back(jet);
   }
}


/**
 * @brief store as _splitting the most recent shared ancestor between the maximum pt particle and antiparticle in the event
 *  
 * this is one definition of the Pythia-level splitting history, which from two particles selected in the final state (in this case, the maximum pt particle and antiparticle)
 * it defines their splitting to be their most recent common ancestor, subject to the requirements that the ancestor is a gluon and the splittees have the same flavor as the particle
 * and were produced in the shower. Compared to find_splitting_v2 (see below) this has the advantage that it is closer to what is accessible experimentally, but the disadvantage that
 * it fails if the two particles in the final state were in fact not produced in a common splitting (for example, because there are two splittings in the event).
 */
void EventCCbar::find_splitting() {

  int current_particle, other_particle;
  // status codes that pythia uses to identify incoming particles of various types (hard event, subevent, and initial state radiation)
  vector<int> incoming_status {constants::INCOMING_HARD, constants::INCOMING_SUB, constants::INCOMING_ISR};

  // current particle and other particle are the two highest-pT tagged particles
  current_particle = _maxpt_tagged_particle.user_info<ExtraInfo>().global_index();
  other_particle = _maxpt_tagged_antiparticle.user_info<ExtraInfo>().global_index();

  // fix one of the particles (other particle) and trace back through the mothers of the other (current particle) until you find a common ancestor, or the beginning of the event
  while ( !_py_event[ other_particle ].isAncestor( current_particle ) && !is_contained( _py_event[ current_particle ].status(), incoming_status ) ) {
    // set current particle to its mother and continue
    current_particle = _py_event[ current_particle ].mother1();
    
  }

  // if the loop has terminated because current particle is an ancestor of other particle 
  // (and its also an ancestor of the other tagged particle, by contruction since we went through its mothers)
  // then the current particle is the originator of the most recent common splitting
  _splitting._in = as_pseudojet( _py_event[ current_particle ] );
  _splitting._out1 = as_pseudojet( _py_event[ _py_event[ current_particle ].daughter1() ] );
  _splitting._out2 = as_pseudojet( _py_event[ _py_event[ current_particle ].daughter2() ] );

  _splitting._is_valid = true;
  // however, some things can go wrong, where we got to the beginning of the event or got some weird splitting, so flag that
  // mother is not a gluon, not a valid splitting
  if ( _splitting._in.user_info<ExtraInfo>().pdg_id()!=constants::GLUON ) _splitting._is_valid = false;
  // daughters are not the identified particle type, not a valid splitting
  if ( abs(_splitting._out1.user_info<ExtraInfo>().pdg_id())!=abs(_particle_ids[0]) || abs(_splitting._out2.user_info<ExtraInfo>().pdg_id())!=abs(_particle_ids[0]) ) _splitting._is_valid = false;
  // daughters were not produced in the shower, not a valid splitting
  if ( _py_event[ _splitting._out1.user_info<ExtraInfo>().global_index() ].status()!=constants::SHOWER || _py_event[ _splitting._out2.user_info<ExtraInfo>().global_index() ].status()==constants::SHOWER ) _splitting._is_valid = false;

  // for later convenience we want _out2 to be the lower pT of the two outgoing particles, so swap them if its not
  if (_splitting._out2.perp() > _splitting._out1.perp() ) swap( _splitting._out1, _splitting._out2 );

  // calculate and set parameters of the splitting from the in and out particles
  _splitting.set_values();

}


/**
 * @brief store as _splitting the highest-pt gluon->ccbar splitting for which the gluon is within jetR of the jet axis
 * 
 * this is one definition of the Pythia-level splitting history, goes through the entire pythia event and find the highest-pt g->ccbar splitting
 * that is within deltaR=jetR of the jet axis. This has the disadvantage compared to find_splitting (above) that it circumvents subtelties in finding
 * a final-state ccbar pair that were produced in the same splitting, but it is perhaps a more correct pythia baseline since it is not affected by these 
 * subtelties.
 * @param jet: the jet under consideration; not that only information about the jet axis is used, for determining if the gluon is inside the jet
 * @param jetR: the allowable distance of the gluon from the jet axis. For most purposes, this should probably be the same as the jet radius
 * @return true if a splitting was found, false otherwise
 */
bool EventCCbar::find_splitting_v2(PseudoJet jet, double jetR) {

  int n_splittings = 0;
  double max_mother_pt = 0.0;

  // loop over the pythia event
  for (Particle& p: _py_event) {

    if (p.status()!=constants::SHOWER) continue; // only consider particles produced through showering
    if (p.id()!=constants::CHARM) continue; // only consider charms

    Particle& mother = _py_event[ p.mother1() ];
    if ( mother.id()!=constants::GLUON ) continue; // only consider particles whose mother is a gluon
    fastjet::PseudoJet g = as_pseudojet( mother );
    if (g.delta_R(jet) < jetR) { // require that the gluon be within 1 unit of the jet axis
      n_splittings++;
      if (mother.pT() > max_mother_pt) { // set the splitting to be the one whose mother has largest pT

        _splitting._in = as_pseudojet( mother );
        _splitting._out1 = as_pseudojet( _py_event[ mother.daughter1()] );
        _splitting._out2 = as_pseudojet( _py_event[ mother.daughter2() ] );
        max_mother_pt = mother.pT();
      }
    }
  }
  if (_splitting._out1.perp() < _splitting._out2.perp()) swap( _splitting._out1, _splitting._out2);
  _splitting.set_values();
  return n_splittings!=0;
}


/**
 * @brief function to calculate the level of a splitting in the pythia shower and set _splitting._level
 * 
 * the level here is defined as the number of splittings with momentum sharing fraction z>0.1 that happened earlier in the shower than _splitting
 */
void EventCCbar::calculate_splitting_level() {

  // status codes that pythia uses to identify incoming particles of various types (hard event, subevent, and initial state radiation)
  vector<int> outgoing_status {constants::OUTGOING_HARD, constants::OUTGOING_SUB, constants::OUTGOING_ISR};

  int current_particle = _py_event[ _splitting._in.user_info<ExtraInfo>().global_index() ].iTopCopy();
  int level = 1, mother_index; // note that this starts with the mother of the ccbar, so level should start at 1 rather than 0
  double ptratio=0.;
  // fix one of the particles (other particle) and trace back through the mothers of the other (current particle) until you find a common ancestor, or the beginning of the event
  while ( !is_contained( _py_event[ current_particle ].status(), outgoing_status ) ) {
    // set current particle to its mother and continue
    mother_index = _py_event[ current_particle ].mother1();

    // increment the level counter only if the splitting involved a momentum sharing of at least 0.1
    ptratio = _py_event[current_particle].pT() / _py_event[mother_index].pT();
    if ( ptratio>0.1 && ptratio<0.9 ) {
      level++;
    }
    // set the current particle to the "top" mother (iTopCopy ignores momentum reshuffling on the same particle)
    current_particle = _py_event[ mother_index ].iTopCopy();
  }
  _splitting._level = level;
}

/**
 * @brief function to iteratively recluster a jet following the maximum pt tagged particle and antiparticle in the reclustering
 * 
 * reclusters a jet iteratively into subjets until the first occasion where the maximum pt particle and antiparticle are split into separate subjets.
 * According to the reclustering logic, this is the reclustered splitting. Set values of _recl_splitting.
 */
void EventCCbar::do_iterative_reclustering(PseudoJet jet) { 

  assert( (contains(jet, _maxpt_tagged_particle) && contains(jet, _maxpt_tagged_antiparticle)) && "jet is expected to contain both the highest pt particle and antiparticle" );

  vector<PseudoJet> jet_constituents = jet.constituents();
  ClusterSequence reclustering(jet_constituents, _jet_def_recl);
  PseudoJet reclustered_jet = (reclustering.inclusive_jets(0.0)).at(0);  // get the reclustered jet
    
  PseudoJet subjet_1, subjet_2;
  bool jet1_contains_particle, jet1_contains_antiparticle, jet2_contains_particle, jet2_contains_antiparticle;
  
  int level=0;
  _recl_splitting._is_primary=true;
  while (reclustered_jet.has_parents( subjet_1, subjet_2)) {

    jet1_contains_particle = contains(subjet_1, _maxpt_tagged_particle);
    jet1_contains_antiparticle = contains(subjet_1, _maxpt_tagged_antiparticle);
    jet2_contains_particle = contains(subjet_2, _maxpt_tagged_particle);
    jet2_contains_antiparticle = contains(subjet_2, _maxpt_tagged_antiparticle);

    // increment the level only if the splitting had z>0.1
    if ( (min(subjet_1.perp(), subjet_2.perp()) / reclustered_jet.perp())>0.1 ) level++; 

    // if each subjet contains exactly one of the particle/ antiparticle pair, and between them they have both, return this splitting
    if ( jet1_contains_particle && jet2_contains_antiparticle ) break;
    if ( jet1_contains_antiparticle && jet2_contains_particle ) break;

    // if the subjet that contains both particles is the softer one at ANY splitting in the reclustering, the ccbar pair is not in the primary Lund plane
    if ( jet1_contains_particle && jet1_contains_antiparticle ) { // subjet 1 contains both charm and anticharm
      if (subjet_1.perp() < subjet_2.perp()) _recl_splitting._is_primary = false;
      // if (subjet_2.perp() > 0.2*reclustered_jet.perp()) _recl_splitting._is_primary = false; // if the jet without the ccbar is too hard, reject
      reclustered_jet = subjet_1; // reiterate with the subjet that contains both particles 
    }
    // same as above except for the case that both charms are in subjet 2
    else if ( jet2_contains_particle && jet2_contains_antiparticle ) {
      if (subjet_2.perp() < subjet_1.perp()) _recl_splitting._is_primary = false;
      // if (subjet_1.perp() > 0.2*reclustered_jet.perp()) _recl_splitting._is_primary = false; // if the jet without the ccbar is too hard, reject
      reclustered_jet = subjet_2;
    }                                           
  }

  // now that each subjet contains on of the particles in tagged_list, define subjet 1 to be the higher pt one
  if ( subjet_2.perp()>subjet_1.perp() ) swap( subjet_1, subjet_2 );

  _recl_splitting._level = level;
  _recl_splitting._in = reclustered_jet;
  _recl_splitting._out1 = subjet_1;
  _recl_splitting._out2 = subjet_2;
  _recl_splitting.set_values();
}


// Splitting EventCCbar::find_hardest_splitting(fastjet::PseudoJet jet) {

//   vector<int> outgoing_status {constants::OUTGOING_HARD, constants::OUTGOING_SUB, constants::OUTGOING_ISR}; // status codes classified by pythia as from incoming particles in the event

//   // find the incoming particle in the event that is an ancestor of at least 90% of the particles in the jet
//   vector<int> initial_outgoing_particles;
//   for (Particle& p: _py_event) {
//     if (p.pT()<10) continue; // don't consider particles with less than 10GeV as potential "initiators"
//     if (is_contained( p.status(), outgoing_status)) initial_outgoing_particles.push_back(p.index());
//   }

// vector<int> count_outgoing(initial_outgoing_particles.size());

// for (int i=0; i<initial_outgoing_particles.size(); i++) {
//   for (auto c: jet.constituents()) { // loop over jet constituents
//     if ( _py_event[c.user_info<ExtraInfo>().global_index()].isAncestor(initial_outgoing_particles[i])) {
//       count_outgoing[i]+=1;
//     }
//   }
// }
// // the initiator is the particle in the event that is the ancestor of the most particles in the jet
// int index = max_element(count_outgoing.begin(),count_outgoing.end()) - count_outgoing.begin();
// int initiator_index = initial_outgoing_particles[index];

// // iterate through splittings of the initiator, following the hardest branch
// int current_particle = _py_event[initiator_index].iBotCopy();
// int daughter1, daughter2;
// double max_kt = -1e9;
// Splitting max_kt_split;

// while (abs(_py_event[current_particle].status())<60) { // before the status code where we are preparing for hadronization

//   daughter1 = _py_event[current_particle].daughter1();
//   daughter2 = _py_event[current_particle].daughter2();
//   if ( _py_event[daughter1].pT() < _py_event[daughter2].pT() ) swap(daughter1, daughter2); // make sure daughter1 is the higher-pt one
//   double kt = get_kt( as_pseudojet(_py_event[daughter1]), as_pseudojet(_py_event[daughter2]) ); // get kt of the current splitting
//   if (kt > max_kt) {
//     max_kt = kt;
//     max_kt_split._in = _py_event[current_particle];
//     max_kt_split._out1 = _py_event[daughter1];
//     max_kt_split._out2 = _py_event[daughter2];
//     max_kt_split._in_index = _py_event[current_particle].index();
//     max_kt_split._out1_index = _py_event[daughter1].index();
//     max_kt_split._out2_index = _py_event[daughter2].index();
//   }
//   current_particle = _py_event[daughter1].iBotCopy(); // set the current particle to be the "most bottom" copy of the daughter (after recoilers)
// }
// return max_kt_split;
// }

