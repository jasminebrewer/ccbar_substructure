#include "histograms.hh"
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace std;

/**
 * @brief function to fill a histogram at location value x with value y
 * 
 * @param xvalue: xvalue from which to determine which histogram bin the addition should go into
 * @param yvalue: the value to add to that histogram bin
*/
void Histogram::fill(double xvalue, double yvalue) {
  
  _values[ getbin(xvalue) ] += yvalue;
}


/**
 * @brief function to add a new histogram to the current histogram
 * 
 * @param new_histogram: histogram whose values to add to the values in the current histogram. Must be the same size as the current histogram.
*/
void Histogram::add(vector<double> new_histogram) {

  assert( new_histogram.size()==_values.size() && "Error: trying to add histograms of different sizes.");
  std::transform(_values.begin(), _values.end(), new_histogram.begin(),
		 _values.begin(), std::plus<double>());
}


/**
 * @brief function to add a new histogram to the current histogram with a multiplicative weight
 * 
 * @param new_histogram: histogram whose values to add to the values in the current histogram. Must be the same size as the current histogram.
 * @param weight: multiplicative weight to apply to new_histogram before adding it to the current histogram
*/
void Histogram::add(vector<double> new_histogram, double weight) {

  assert( new_histogram.size()==_values.size() && "Error: trying to add histograms of different sizes.");
  // first multiply the new histogram values by the weight
  std::transform(new_histogram.begin(), new_histogram.end(), new_histogram.begin(),
		 std::bind1st(std::multiplies<double>(), weight));
  // then add them to the histogram
  std::transform(_values.begin(), _values.end(), new_histogram.begin(),
		 _values.begin(), std::plus<double>());
}


/**
 * @brief function to multiply a histogram by a constant multiplier
 * 
 * @param multiplier: value to multiply the histogram values by
*/
void Histogram::multiply(double multiplier) {

  std::transform(_values.begin(), _values.end(), _values.begin(),
		 std::bind1st(std::multiplies<double>(), multiplier));
}


/**
 * @brief function to get the appropriate bin for a value in the histogram
 * 
 * @param value: x value to get the bin for
 * @return integer index of the bin (out of bounds are stored in location 0, as with np.digitize)
*/
int Histogram::getbin( double value ){

  double set_point, bin_size;
  set_point = _bins[1]; // minimum values are in place one (so out of bounds can be stored in 0) (same convention as np.digitize)
  bin_size = _bins[2]-_bins[1];
  
  if (value < _bins[1] || isinf(value) || isnan(value) ) return 0;
  else if (value>=_bins[_bins.size()-2]) return _bins.size()-2;// since bins is 1 element longer than the eec vectors
  else return int( floor( (value-set_point)/bin_size ) )+1;
}


/**
 * @brief function to write a histogram to _outfile as bin, value pairs
*/
void Histogram::write() {

  for (long unsigned int i=0; i<_values.size(); i++) {
    _outfile << _bins[i] << ", " << _values[i] << endl;
  }
}


/**
 * @brief function to set all values in a histogram to zero
*/
void Histogram::set_to_zero() {

  std::fill(_values.begin(),_values.end(),0.0);
}
