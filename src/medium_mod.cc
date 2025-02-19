#include "complex_Ei.hh"
#include "constants.hh"
#include "medium_mod.hh"
#include "global_event_analysis.hh"
#include <complex>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <fstream>

using namespace std;
const std::complex<double> i(0, 1);


/* three helper functions to compute quantities used by the integrand function */
double c(double z) {
  return 1. - (constants::CA/constants::CF)*z*(1.-z);
}

std::complex<double> compute_omega(double eg, double z) { return sqrt(c(z) / (eg*z*(1.-z)*i) );}

double compute_mu(double eg, double z) { return 0.5*eg*z*(1.-z);}


/**
 * @brief function to compute the integrand whose integral over x\in[0,1] gives the medium modification of the splitting function
 * 
 * @param x: the integration variable
 * @param p: the struct splitting_params that specifies physical parameters for the integrand
 * @return the (double) value of the integrand at x
 */
double integrand(double x, void * p) {

  // retrieve the parameters for the splitting
  struct splitting_params * params = (struct splitting_params *)p;
  double z = (params->z);
  double eg = 2.*(params->Eg) / ((params->L)*(params->qL));
  double mc2t = (params->mc2) / (params->qL);
  double pt2t = (params->pt2) / (params->qL);
  double qL = (params->qL);

  // compute omega and mu values for use in calculating the integrand
  std::complex<double> omega = compute_omega(eg,z);
  double mu = compute_mu(eg,z);

  std::complex<double> exparg1 = -i*0.5*pt2t*tan(omega*x)/(mu*omega);
  std::complex<double> exparg2 = (1.+i*0.5*tan(omega*x)*c(z)*(1.-x)/(mu*omega));
  // integrand involves difference between to Exponential Integral Ei functions, whose arguments in practice are very similar
  // (exparg2 \approx 1 + epsilon*c, where epsilon is small and c is complex)
  // uses custom implementation of the exponential integral for complex arguments which is accurate over wide range of argument values
  std::complex<double> exp_diff = Ei(exparg1)-Ei(exparg1/exparg2);

  std::complex<double> I4 = (1./(qL*eg*c(z))) * exp(-i*0.5*x*mc2t/mu) * (-i*omega/sin(omega*x)) * ( (mc2t-2.*i*mu*omega*(z*z+(1.-z)*(1.-z))/(sin(omega*x)))*exp_diff
											       -(2.*i*mu*omega*(z*z+(1.-z)*(1.-z))/(sin(omega*x)))*(-exp(-i*0.5*pt2t*tan(omega*x)/(mu*omega))
																		    + exp(exparg1/exparg2)/exparg2) );

  std::complex<double> I5 = (1./(qL*eg)) * exp(-i*0.5*mc2t*x/mu) * (1./(i*(mc2t+pt2t)*cos(omega*x))) * exp(-i*0.5*pt2t*tan(omega*x)/(mu*omega)) * (mc2t + (pt2t/cos(omega*x))*(z*z + (1.-z)*(1.-z)));

  // the integrand is the real part of I4+I5
  return (I4 + I5).real();
}


/**
 * @brief function to compute the medium modified splitting by performing the integration over the integrand specified above
 * 
 * @param p: the struct splitting_params that specifies physical parameters for the integrand
 * @param gauss_integrator: boolean for whether you want gauss_kronrod integration (true) or tanh_sinh integration (false). The recommended value is true
 * @return the (double) value of the medium modified splitting function
 */
double med_splitting(void * p, bool gauss_integrator) {

  double result, error, L1;
  struct splitting_params params = * (struct splitting_params *)p;

  auto f = [&params](double x) {
    return integrand(x, &params);
  };
  double termination = 1.e-10;//std::sqrt(std::numeric_limits<double>::epsilon());

  // the gauss_kronrod integration works typically better than the tanh_sinh, which fails epically on (relatively rare) occasion
  if (!gauss_integrator) {
    size_t max_halvings = 200;
    boost::math::quadrature::tanh_sinh<double> integrator(max_halvings);
    result = integrator.integrate(f, 0.0, 1.0, termination, &error, &L1);
  }
  else {
    // result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, 0.0, 1.0, 10, termination, &error);
    result = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, 0.0, 1.0, 15, termination, &error);
  }
  return result;
}


/**
 * @brief function to compute the integrand whose integral over x\in[0,1] gives the medium modification of the splitting function
 * 
 * @param x: the integration variable
 * @param p: the struct splitting_params that specifies physical parameters for the integrand
 * @return the (double) value of the integrand at x
 */
double vac_splitting(void * p) {

  // retrieve the parameters for the splitting
  struct splitting_params * params = (struct splitting_params *)p;
  double z = (params->z);
  double mc2 = (params->mc2);
  double pt2 = (params->pt2);

  // return the vacuum splitting function for those parameters
  return 0.5* ((z*z*(1.-z)*(1.-z))/pow(mc2+pt2,2.)) * (((mc2+pt2)*(z*z+(1.-z)*(1.-z)))/(z*(1.-z)) + 2*mc2);
}

/**
 * @brief compute the medium weight from the medium modified splitting, which is 1 + med_splitting / vac_splitting
 * 
 * @param p: the struct splitting_params that specifies physical parameters for the integrand
 * @param gauss_integrator: whether to do gauss_kronrod integration (see med_splitting for further details)
 * @return weight of the medium event
 */
double compute_medium_weight(void * p, bool gauss_integrator) {

  return 1. + med_splitting(p, gauss_integrator) / vac_splitting(p);

}

/**
 * @brief auxiliary function that can be used to output a grid of the medium weights for a grid in k^2 and z. 
 * 
 * Used for making 2d plots of the splitting function and cross-checking with Mathematica implementation
 */
void make_2d_plot() {

  ofstream outfile;
  outfile.open("Pmed_over_Pvac.txt");

  splitting_params params;
  double invGeVtofm = 0.1973; // conversion for length units in GeV^{-1} to fermi
  params.L = 4.0 / invGeVtofm;
  params.qL = 2.0;
  params.mc2 = pow(1.27, 2.);
  params.Eg = 100.;

  for (int iz=1; iz<51; iz++) {

    for (int ik=1; ik<51; ik++) {

      params.z = 0.5*(iz/50.);
      params.pt2 = 0.2*ik;
      outfile << params.z << ", " << params.pt2 << ", " << compute_medium_weight(&params, true) << endl;
    }
  }
}
