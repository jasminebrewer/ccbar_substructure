#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/gamma.hpp>
#include "constants.hh"
#include "ccbar_analysis.hh"
//#include <unordered_set>

using namespace fastjet;
using namespace std;

struct energy_loss_params {double qL; double L; double n; double T; double alpha_med; double jetR;}; 

bool is_quark( int pid ) {
    if ( abs(pid) >=1 && abs(pid) <=6 ) return true;
    else return false;
}
bool is_quark( PseudoJet particle ) { return is_quark( particle.user_info<ExtraInfo>().pdg_id() ); }

bool is_gluon( int pid ) {

    if ( pid==constants::GLUON ) return true;
    else return false;
}
bool is_gluon( PseudoJet particle ) { return is_gluon( particle.user_info<ExtraInfo>().pdg_id() ); }

int reduce_flavor( int flavor1, int flavor2 ) {

    // 2 gluons make a gluon, a quark-antiquark pair makes a gluon, a quark and a gluon makes a quark
    // no defined addition operatior for two quarks, two antiquarks, or for quarks of different flavors
    
    if ( flavor1==-flavor2 ) return constants::GLUON; // quark-antiquark pair makes a gluon
    if ( is_gluon(flavor1) && is_gluon(flavor2) ) return constants::GLUON; // two gluons makes a gluon
    if ( is_gluon(flavor1) && is_quark(flavor2) ) return flavor2;
    if ( is_quark(flavor1) && is_gluon(flavor2) ) return flavor1;
    else return 0; // error code; need to keep both particles and reduce the flavor later
}


PseudoJet scaleMomentum(PseudoJet particle, double scale) {
  double new_p2 = pow(scale, 2.) * (particle.px()*particle.px() + particle.py()*particle.py() + particle.pz()*particle.pz());
  double new_energy = sqrt( new_p2 + particle.m()*particle.m() );
  PseudoJet new_particle = PseudoJet(scale*particle.px(), scale*particle.py(), scale*particle.pz(), new_energy);
  new_particle.set_user_info(new ExtraInfo(particle.user_info<ExtraInfo>().pdg_id(), particle.user_info<ExtraInfo>().global_index()));
  return new_particle;
}


std::vector<int> remove_pairs(const std::vector<int>& input) {

    vector<int> unpaired_values;

    for (int num : input) {
        auto antiflavor_it = std::find(unpaired_values.begin(), unpaired_values.end(), -num);
        if (antiflavor_it != unpaired_values.end()) {
            // if the opposite value exists, remove it from the list
            unpaired_values.erase( antiflavor_it );
        } else {
            // otherwise, add the current value to the set
            unpaired_values.push_back(num);
        }
    }

    // Convert the unordered_set back to a vector
    return std::vector<int>(unpaired_values.begin(), unpaired_values.end());
}

int compute_flavor(vector<int> flavors) {

    // any number of gluons can be reduced to a single gluon, since two gluons makes a gluon
    // therefore you only need to count unpaired quarks

    // remove all the gluons from the list
    flavors.erase(std::remove(flavors.begin(), flavors.end(), constants::GLUON), flavors.end());

    // remove all matched quark/antiquark pairs
    vector<int> new_flavors = remove_pairs(flavors);

    // if (new_flavors.size()>1) {
    //     cout << "WARNING: more than one flavor. Keeping the first..." << endl;
    //     cout << "flavors: ";
    //     for (auto f: new_flavors) cout << f << ", ";
    //     cout << endl;
    // }

    if (new_flavors.size()==0) return constants::GLUON; // the flavor is gluon since all the quarks are paired
    else return new_flavors[0];
}


double eloss_integrand( double x, void * params) {
    double alpha = *(double *) params;
    double alpha_x2 = alpha*x*x;
    return ( (1.-exp(-x))/pow(x,1.5) )*(exp(-alpha_x2) - alpha_x2*boost::math::tgamma(1e-10,alpha_x2));
};

double Ifun( double alpha, double minval ) {

    double result, error;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    gsl_function F;
    F.function = &eloss_integrand;
    F.params = &alpha;

    gsl_integration_qags (&F, minval, 10.*minval, 0, 1.e-7, 1000,
                        w, &result, &error);

    gsl_integration_workspace_free (w);

    return result;
}

double compute_quenching_weight(int flavor, double pT, void * p) {

    // retrieve energy loss parameters
    struct mediumParams * params = (struct mediumParams *)p;
    double qhatL = params->qL;
    double jetR = params->jetR;
    double n = params->n;
    double T = params->T;
    double alpha_med = params->alpha_med;

    double L = params->L;
    double omega_c = qhatL* L / 2.;

    double omega_s = pow( alpha_med*constants::Nc / M_PI, 2. ) * omega_c;
    double alpha = pow( 2.*pT*jetR / (sqrt(qhatL) * n), 2. );
    
    double Q_minijet, Q_pert;

    double Ci;
    if ( is_quark(flavor) ) Ci=constants::CF;
    else if ( is_gluon(flavor) ) Ci=constants::CA;

    Q_minijet = exp( (-2.*alpha_med*Ci / M_PI) * (sqrt(2.*omega_c/T)*(1.-exp(-n*T/pT)) - sqrt(2.*omega_c/omega_s)*(1.-exp(-n*omega_s/pT))
        + sqrt(2.*M_PI*omega_c*n/pT)*( erf(sqrt(omega_s*n/pT)) - erf(sqrt(n*T/pT)))) );
 
    Q_pert = exp( -(alpha_med*Ci/M_PI) * sqrt(2.*omega_c*n/pT) * Ifun( alpha, n*omega_s/pT ) );

    return Q_minijet * Q_pert;
}

vector<PseudoJet> compute_jet_modification( vector<PseudoJet> jet_constituents, void * p ) {

    // retrieve energy loss parameters
    struct mediumParams * params = (struct mediumParams *)p;
    double qhatL = params->qL;
    double L = params->L;
    double theta_c = 2. / sqrt( qhatL*L*L );
    double jetR = params->jetR;
    string resolution_mode = params->resolutionMode;

    double minijetR;
    if (resolution_mode=="jet") minijetR=jetR; // resolution length of the entire jet, meaning the jet loses energy as a whole object
    else if (resolution_mode=="thetac") minijetR=theta_c; // resolution length of theta_c
    else if (resolution_mode=="zero") minijetR=1e-10; // zero resolution length, so each particle is modified individually
        
    // cluster particles into groups depending on which ones are closer to each other than the medium resolution length
    JetDefinition mini_jet_def(kt_algorithm, minijetR);
    ClusterSequence cs(jet_constituents, mini_jet_def);
    vector<PseudoJet> mini_jets = (cs.inclusive_jets(0.0));  // get all clusters with radius theta_c

    vector<PseudoJet> modified_jet;

    double quenching_weight = 0.;
    for (auto cluster: mini_jets) {
        vector<int> flavors;
        int final_flavor;

        for (auto c: cluster.constituents()) flavors.push_back(c.user_info<ExtraInfo>().pdg_id());
        final_flavor = compute_flavor(flavors);

        // compute the quenching weight of the cluster
        quenching_weight = compute_quenching_weight(final_flavor, cluster.perp(), p);

        for (auto c: cluster.constituents()) {
            modified_jet.push_back( scaleMomentum(c, quenching_weight) );
        }
    }

    return modified_jet;

}