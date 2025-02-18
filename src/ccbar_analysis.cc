#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "ccbar_analysis.hh"
#include "constants.hh"
#include "histograms.hh"
#include "splitting.hh"
#include "jet_energy_loss.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/FlavorCone.hh"


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

PseudoJet as_pseudojet( vector<PseudoJet> particles ) {

  JetDefinition jet_def(fastjet::antikt_algorithm, 5.0); // very large jet radius ensures all particles are within the jet radius
  ClusterSequence cs(particles, jet_def);
  // store all jets in the event above the min pT cut and sort them by pT
  vector<PseudoJet> all_jets = sorted_by_pt( cs.inclusive_jets(0.0) );
  
  if (all_jets.size()!=1) cout << "ERROR in as_pseudojet: incorrect number of jets." << endl;
  return all_jets[0];
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
 * @brief helper function to return an empty PseudoJet object, with ExtraInfo defined
 */
PseudoJet emptyPseudoJet() {
  PseudoJet empty = PseudoJet(0.,0.,0.,0.);
  empty.set_user_info(new ExtraInfo(0, -1));
  return empty;
}

/**
 * @brief helper function to associate a particle with one of a list of other particles, based on which one has the minimum deltaR distance
 *  
 * primary use case is to associate a parton with its closest hadron
 * 
 * @param particle PseudoJet object referring to one particle
 * @param tagged_particles vector of PseudoJets referring to a number of particles
 * @return the particle in tagged_particles that is closest in deltaR to particle
 */
PseudoJet associate(PseudoJet particle, vector<PseudoJet> tagged_particles) {

  double min_dR=1.e10;
  PseudoJet associated_particle;
  
  for (auto tp: tagged_particles) {
    // both particles have to have the same sign of their pdg (i.e. either both are particles, or both are antiparticles)
    if ( signbit(tp.user_info<ExtraInfo>().pdg_id())!=signbit(particle.user_info<ExtraInfo>().pdg_id()) ) continue;
    // then keep the one with minimum deltaR
    if ( tp.delta_R(particle) < min_dR ) {
      associated_particle = tp;
      min_dR = tp.delta_R(particle);
    }
  }
  if ( min_dR>0.4 ) associated_particle = emptyPseudoJet(); // association failed!
  return associated_particle;
}

/**
 * @brief function to trace a non-final state pythia particle to the corresponding final state particle
 *  
 * the particle is assumed to have a unique daughter with the same particle id which allows that particle to be traced to the final state
 * (namely, this function should be used for particles undergoing q->qg splittings, and the function follows the q)
 * 
 * @param particle PseudoJet object for a particular particle to trace to the final state (should be a parton)
 * @return PseudoJet object associated with the final state version of that particle, which is either a parton or a hadron depending on 
 * whether the event is parton- or hadron-level
 */
PseudoJet EventCCbar::follow_to_final_state( PseudoJet particle) {

  int current_index = particle.user_info<ExtraInfo>().global_index(); // event index of the particle
  int pid = particle.user_info<ExtraInfo>().pdg_id(); // its particle id
  Particle& current_particle = _py_event[ current_index ]; // pythia particle associated with that particle

  // while the current particle is not final
  while ( (_is_parton_level && !current_particle.isFinal()) || (!_is_parton_level && !current_particle.isFinalPartonLevel()) ) {
 
    // set the new current particle to which daughter has the same particle id as the current particle
    if ( _py_event[ current_particle.daughter1() ].id()==pid) current_index = current_particle.daughter1();
    else if ( _py_event[ current_particle.daughter2() ].id()==pid) current_index = current_particle.daughter2();

    current_particle = _py_event[ current_index ];
  }

  // if the event is hadron-level, associate the particle with one of the tagged hadrons in the event
  if (!_is_parton_level) return associate( as_pseudojet(_py_event[current_index]), _tagged_particles );
  // else convert the final particle to a PseudoJet object and return
  else return as_pseudojet( _py_event[current_index] );
}


void EventCCbar::add_daughter( PseudoJet particle, vector<PseudoJet>& daughters) {

  int current_index = particle.user_info<ExtraInfo>().global_index(); // event index of the particle
  int pid = particle.user_info<ExtraInfo>().pdg_id(); // its particle id
  Particle& current_particle = _py_event[ current_index ]; // pythia particle associated with that particle

  for (auto p: _final_particles) {
    int cid = p.user_info<ExtraInfo>().global_index();
    Particle& pp = _py_event[ cid ];
    // cout << "outside and is final? " << pp.isFinal() << endl;
    if ( pp.isFinal() && pp.isAncestor(current_index) ) {
      daughters.push_back(p);
    }
  }

  // for (auto p: current_particle.daughterListRecursive()) cout << _py_event[ p ].index() << ", ";

  // if ( (_is_parton_level && !current_particle.isFinal()) || (!_is_parton_level && !current_particle.isFinalPartonLevel()) ) {

  //   PseudoJet daughter1 = as_pseudojet( _py_event[current_particle.daughter1()] );
  //   PseudoJet daughter2 = as_pseudojet( _py_event[current_particle.daughter2()] );
  //   add_daughter( daughter1, daughters );
  //   add_daughter( daughter2, daughters );
  // }
  // if ( (_is_parton_level && current_particle.isFinal()) || (!_is_parton_level && current_particle.isFinalPartonLevel())) {
  //   daughters.push_back(particle);
  // }
}
PseudoJet EventCCbar::all_daughters_until_final_state( PseudoJet particle) {

  vector<PseudoJet> daughters;
  add_daughter( particle, daughters );
  return as_pseudojet( daughters );
}


/**
 * @brief read in a pythia event
 *  
 * from the pythia event, store final state particles passing the track selection cuts, and store global indices of particles of particular interest to the analysis, as class members
 */
void EventCCbar::read_event() {

  // vector<int> cids = {4,-4};
  // vector<int> bids = {5,-5};

  // loop over all particles in the event
  for (Particle& p: _py_event) {

    // perform cuts to remove particles that would have gone out of the detector
    if (p.pT() < _track_cuts.trackPtMin)  continue;
    if (abs(p.eta()) > _track_cuts.trackEtaCut) continue;
    
    // if hadronization is on, store tagged particles at parton-level to associate them later with the tagged hadrons
    if (!_is_parton_level && p.isFinalPartonLevel()) {
      if (is_contained(p.id(), _parton_ids) && p.pT()>_track_cuts.HFPtMin) {
        PseudoJet parton( p.px(), p.py(), p.pz(), p.e() );
    
        // store the particle id and global index as part of about the particle
        parton.set_user_info(new ExtraInfo(p.id(), p.index()));
        _tagged_partons.push_back( parton );
      }
    }

    if (!p.isFinal()) continue; // otherwise, consider only final states

    
    PseudoJet particle( p.px(), p.py(), p.pz(), p.e() );
    
    // store the particle id and global index as part of about the particle
    particle.set_user_info(new ExtraInfo(p.id(), p.index()));
    
    _final_particles.push_back(particle);

    // if (is_contained(p.id(), cids)) _tagged_cquarks.push_back(particle);
    // if (is_contained(p.id(), bids)) _tagged_bquarks.push_back(particle);

    if (_is_inclusive) continue; // no extra particle info to save, so we're done

    // add a particle to tagged_particles if it has a particle id matching particle_ids, and if its pT is bigger than the low pT cut for heavy flavor 
    if (is_contained(p.id(), _particle_ids) && particle.perp()>_track_cuts.HFPtMin ) {
      _tagged_particles.push_back(particle);
    }
  }
  // sort the tagged particles by pT
  _tagged_particles = sorted_by_pt(_tagged_particles); 
  _has_pair = get_pair(); // set values for maxpt_tagged_particle and maxpt_tagged_antiparticle, and has_pair

}


/**
 * @brief cluster all final-state particles in the event into jets
 *  
 * cluster jets according to the provided jet definition, reject jets that don't satisfy the jet selection criteria, and store up to two highest pt jets as class members
 */void EventCCbar::cluster_jets() {

   // cluster final particles in the event according to the jet definition jet_def and save the associated cluster sequence for later
   vector<PseudoJet> all_jets;

   ClusterSequence cs_umod(_final_particles, _jet_def);
   // store all jets in the event above the min pT cut and sort them by pT
   all_jets = sorted_by_pt( cs_umod.inclusive_jets(_track_cuts.JetPtMin) );

   // if (all_jets.size() > 0) cout << "unmodified jets: ";
   for (auto jet: all_jets) {

     // reject jets that either don't have constituents or fall outside of the acceptance cuts of the analysis
     bool pass_jet_cuts = (jet.has_constituents() && jet.eta()>_track_cuts.JetEtaMin && jet.eta()<_track_cuts.JetEtaMax && jet.pt()>_track_cuts.JetPtMin && jet.pt()<_track_cuts.JetPtMax );
     if (!pass_jet_cuts) continue;

     // reject jets that should have a c-cbar pair, but don't
     // if (!_is_inclusive && !get_pair(jet)) continue;

     if ( _unmodified_jets.size() < 2 ) _unmodified_jets.push_back(jet);
   }

   // to save computational time, only compute energy loss if at least some (unmodified) jets pass the cuts and satisfy the event selection criteria for the analysis
   bool is_candidate_event = ( _unmodified_jets.size()>0 && (_has_pair || _is_inclusive) );

   if (_do_energy_loss && is_candidate_event) {
    _final_particles = compute_jet_modification( _final_particles, &_medium_params);

    // update tagged particles and maxpt particle and antiparticle
    _tagged_particles.clear();
    for (auto p: _final_particles) {
      // add a particle to tagged_particles if it has a particle id matching particle_ids, and if its pT is bigger than the low pT cut for heavy flavor 
      if (is_contained(p.user_info<ExtraInfo>().pdg_id(), _particle_ids) && p.perp()>_track_cuts.HFPtMin ) {
        _tagged_particles.push_back(p);
      }
    }
    // sort the tagged particles by pT
    _tagged_particles = sorted_by_pt(_tagged_particles);
   }
   _has_pair = get_pair(); // reset _maxpt_tagged_particle and _maxpt_tagged_antiparticle

    ClusterSequence cs(_final_particles, _jet_def);
    _cluster_seq = cs;
    // store all jets in the event above the min pT cut and sort them by pT
    all_jets = sorted_by_pt( _cluster_seq.inclusive_jets(_track_cuts.JetPtMin) );

    // if (all_jets.size() > 0) cout << "modified jets: ";
    for (auto jet: all_jets) {

      // reject jets that either don't have constituents or fall outside of the acceptance cuts of the analysis
      bool pass_jet_cuts = (jet.has_constituents() && jet.eta()>_track_cuts.JetEtaMin && jet.eta()<_track_cuts.JetEtaMax && jet.pt()>_track_cuts.JetPtMin && jet.pt()<_track_cuts.JetPtMax );
      if (!pass_jet_cuts) continue;

      // reject jets that should have a c-cbar pair, but don't
      // if (!_is_inclusive && !get_pair(jet)) continue;

      // add at most 2 highest pt jets from the event satisfying the cuts to the event's jets
      if ( _jets.size() < 2 ) _jets.push_back(jet);
   }
}




std::pair<int,int> EventCCbar::find_typical_initiator(PseudoJet jet) {

  unordered_map<int, float> values;
  int current_particle;
  for (auto p: sorted_by_pt(jet.constituents())) {
    current_particle = p.user_info<ExtraInfo>().global_index();
    PseudoJet initiator = find_initiator(current_particle);
    values[ initiator.user_info<ExtraInfo>().global_index() ] += p.perp();
  }

  auto maxPair = *std::max_element(
    values.begin(), values.end(),
    [](const auto& a, const auto& b) { return a.second < b.second; }
  );

  return {_py_event[maxPair.first].id(), _py_event[maxPair.first].status()};

}


//   // status codes that pythia uses to identify incoming particles of various types (hard event, subevent, and initial state radiation)
//   vector<int> incoming_status {constants::INCOMING_HARD, constants::INCOMING_SUB, constants::INCOMING_ISR};

//   int current_particle, mother_particle;
//   // current particle and other particle are the two highest-pT tagged particles
//   current_particle = _maxpt_tagged_particle.user_info<ExtraInfo>().global_index();
//   mother_particle = _py_event[ current_particle ].mother1();

//   int i =0;

//   cout << "inside: " << endl;
//   // fix one of the particles (other particle) and trace back through the mothers of the other (current particle) until you find a common ancestor, or the beginning of the event
//   while ( !is_contained( _py_event[ mother_particle ].status(), incoming_status ) && i<20 ) {
//     // set current particle to its mother and continue
//     cout << _py_event[current_particle].status() << ", ";
//     current_particle = _py_event[ current_particle ].mother1();
//     mother_particle = _py_event[ current_particle ].mother1();
//     i++;
//   }
//   cout << endl;
//   if (i==20) cout << "large i! " << _py_event[current_particle].status() << endl;
//   return _py_event[current_particle];
// }

PseudoJet EventCCbar::find_initiator(int current_particle) {

  // status codes that pythia uses to identify incoming particles of various types (hard event, subevent, and initial state radiation)
  vector<int> incoming_status {constants::INCOMING_HARD, constants::INCOMING_SUB, constants::INCOMING_ISR, constants::INCOMING_BEAM_REMNANT};

  int mother_particle;
  // current particle and other particle are the two highest-pT tagged particles
  // current_particle = _maxpt_tagged_particle.user_info<ExtraInfo>().global_index();
  mother_particle = _py_event[ current_particle ].mother1();

  // fix one of the particles (other particle) and trace back through the mothers of the other (current particle) until you find a common ancestor, or the beginning of the event
  while ( !is_contained( _py_event[ mother_particle ].status(), incoming_status ) && _py_event[current_particle].status()!=-11 ) {
    // set current particle to its mother and continue
    current_particle = _py_event[ current_particle ].mother1();
    mother_particle = _py_event[ current_particle ].mother1();
  }

  return as_pseudojet(_py_event[current_particle]);
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
bool EventCCbar::find_splitting_v2(PseudoJet jet, double jetR, bool include_recursive_daughters) {

  int n_splittings = 0;
  double max_mother_pt = 0.0;

  // loop over the pythia event
  for (Particle& p: _py_event) {

    if (p.status()!=constants::SHOWER) continue; // only consider particles produced through showering
    if (!is_contained(p.id(), _parton_ids)) continue; // only consider particles of the type you are trying to tag

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

  // as a test, include all particles from the subsequent radiations in the MC-level definition
  if (include_recursive_daughters) {
    _splitting._in = all_daughters_until_final_state(_splitting._in);
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

bool EventCCbar::get_pair(PseudoJet jet) {

bool has_particle = false;
bool has_antiparticle = false;
double max_particle_pt = 0.0;
double max_antiparticle_pt = 0.0;

for (auto tp: _tagged_particles) {
  if (contains(jet, tp) && tp.user_info<ExtraInfo>().pdg_id()>0 && tp.perp()>max_particle_pt ) {
    _maxpt_tagged_particle = tp;
    max_particle_pt = tp.perp();
    has_particle = true;
  }
  if (contains(jet, tp) && tp.user_info<ExtraInfo>().pdg_id()<0 && tp.perp()>max_antiparticle_pt) {
    _maxpt_tagged_antiparticle = tp;
    max_antiparticle_pt = tp.perp();
    has_antiparticle = true;
  }
}
return (has_particle && has_antiparticle);
}

bool EventCCbar::get_pair() {

bool has_particle = false;
bool has_antiparticle = false;
double max_particle_pt = 0.0;
double max_antiparticle_pt = 0.0;

for (auto tp: _tagged_particles) {
  if (tp.user_info<ExtraInfo>().pdg_id()>0 && tp.perp()>max_particle_pt ) {
    _maxpt_tagged_particle = tp;
    max_particle_pt = tp.perp();
    has_particle = true;
  }
  if (tp.user_info<ExtraInfo>().pdg_id()<0 && tp.perp()>max_antiparticle_pt) {
    _maxpt_tagged_antiparticle = tp;
    max_antiparticle_pt = tp.perp();
    has_antiparticle = true;
  }
}
return (has_particle && has_antiparticle);
}

void record_splittings(PseudoJet reclustered_jet, vector<Splitting>& splitting_list) {

  PseudoJet subjet_1, subjet_2;
  if (reclustered_jet.has_parents( subjet_1, subjet_2 )) {

    Splitting split( reclustered_jet, subjet_1, subjet_2 );
    split.set_values();

    // if (split._virt > 2.*constants::CHARM_MASS) splitting_list.push_back(split);
    if (split._kt > 1.) splitting_list.push_back(split);
    record_splittings(subjet_1, splitting_list);
    record_splittings(subjet_2, splitting_list);
  }
}

Splitting EventCCbar::get_random_splitting(PseudoJet jet) {

  // recluster the jet
  vector<PseudoJet> jet_constituents = jet.constituents();
  ClusterSequence reclustering(jet_constituents, _jet_def_recl);
  PseudoJet reclustered_jet = (reclustering.inclusive_jets(0.0)).at(0);  // get the reclustered jet

  vector<Splitting> all_splittings;
  record_splittings(reclustered_jet, all_splittings);

  if (all_splittings.size()==0) {
    Splitting empty;
    empty.set_values();
    return empty;
  }

  // randomly choose one of the splittings passing the (charm) mass threshold
  int rand_index = rand() % all_splittings.size();
  return all_splittings[ rand_index ];
}


/**
 * @brief function to iteratively recluster a jet following the maximum pt tagged particle and antiparticle in the reclustering
 * 
 * reclusters a jet iteratively into subjets until the first occasion where the maximum pt particle and antiparticle are split into separate subjets.
 * According to the reclustering logic, this is the reclustered splitting. Set values of _recl_splitting.
 */
Splitting EventCCbar::do_iterative_reclustering(PseudoJet jet) { 

  assert( (contains(jet, _maxpt_tagged_particle) && contains(jet, _maxpt_tagged_antiparticle)) && "jet is expected to contain both the highest pt particle and antiparticle" );

  Splitting recl_splitting;

  vector<PseudoJet> jet_constituents = jet.constituents();
  ClusterSequence reclustering(jet_constituents, _jet_def_recl);
  PseudoJet reclustered_jet = (reclustering.inclusive_jets(0.0)).at(0);  // get the reclustered jet
    
  PseudoJet subjet_1, subjet_2;
  bool jet1_contains_particle, jet1_contains_antiparticle, jet2_contains_particle, jet2_contains_antiparticle;
  
  int level=0;
  recl_splitting._is_primary=true;
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
      if (subjet_1.perp() < subjet_2.perp()) recl_splitting._is_primary = false;
      // if (subjet_2.perp() > 0.2*reclustered_jet.perp()) _recl_splitting._is_primary = false; // if the jet without the ccbar is too hard, reject
      reclustered_jet = subjet_1; // reiterate with the subjet that contains both particles 
    }
    // same as above except for the case that both charms are in subjet 2
    else if ( jet2_contains_particle && jet2_contains_antiparticle ) {
      if (subjet_2.perp() < subjet_1.perp()) recl_splitting._is_primary = false;
      // if (subjet_1.perp() > 0.2*reclustered_jet.perp()) _recl_splitting._is_primary = false; // if the jet without the ccbar is too hard, reject
      reclustered_jet = subjet_2;
    }                                           
  }

  // now that each subjet contains on of the particles in tagged_list, define subjet 1 to be the higher pt one
  if ( subjet_2.perp()>subjet_1.perp() ) swap( subjet_1, subjet_2 );

  recl_splitting._level = level;
  recl_splitting._in = reclustered_jet;
  recl_splitting._out1 = subjet_1;
  recl_splitting._out2 = subjet_2;
  recl_splitting.set_values();
  return recl_splitting;
}


Splitting EventCCbar::do_flavor_cone(string mode) {

	Splitting splitting_cc;

  if (mode=="bare") {

  	// for comparison, also include the "splitting" which is just the kinematics of the ccbar pair, without anything fancy
	  splitting_cc._in = _maxpt_tagged_particle+_maxpt_tagged_antiparticle;
	  if (_maxpt_tagged_particle.perp() > _maxpt_tagged_antiparticle.perp() ) {
	    splitting_cc._out1 = _maxpt_tagged_particle;
	    splitting_cc._out2 = _maxpt_tagged_antiparticle;
	  }
	  else {
	    splitting_cc._out1 = _maxpt_tagged_antiparticle;
	    splitting_cc._out2 = _maxpt_tagged_particle;              
	  }
  }
  else if (mode=="standard" || mode=="HFradius") {

    double FC_jetR;
    if (mode=="standard") FC_jetR = 0.5 * _track_cuts.jetR;
    else if (mode=="HFradius") FC_jetR = 0.5 * _maxpt_tagged_particle.delta_R(_maxpt_tagged_antiparticle);

    vector<fastjet::PseudoJet> seeds {_maxpt_tagged_particle, _maxpt_tagged_antiparticle};
    fastjet::contrib::FlavorConePlugin jdf(seeds, FC_jetR);
    fastjet::ClusterSequence jcs(_final_particles, &jdf); 
    vector<fastjet::PseudoJet> FC_jets = sorted_by_pt( jcs.inclusive_jets(0) ); // note that with two seeds, there are always two jets with the FlavorCone method

	  splitting_cc._in = (FC_jets[0] + FC_jets[1]);
	  splitting_cc._out1 = FC_jets[0];
	  splitting_cc._out2 = FC_jets[1];
  }
	splitting_cc.set_values();
  return splitting_cc;
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

