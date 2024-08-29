#ifndef __HISTOGRAMS_HH__
#define __HISTOGRAMS_HH__

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cassert>

using namespace std;

class Histogram {

public:

  Histogram() = default;

  Histogram(double min_bin, double max_bin, double bin_size) {

    // initialize bins according to min, max and bin size provided
    double bin_edge=min_bin;
    while (bin_edge<max_bin) {
      _bins.push_back(bin_edge);
      bin_edge += bin_size;
    }
    if (abs(_bins.back()-max_bin)>1e-10) _bins.push_back(max_bin);
    
    // initialize values with zeros. Note that _values is one element longer than _bins, since we keep a bin both for underflow (leftmost) and for overflow (rightmost)              
    _values.resize(_bins.size()+1);
    std::fill(_values.begin(),_values.end(),0.0);
  }

  Histogram(double min_bin, double max_bin, double bin_size, string filename) {

    // initialize bins according to min, max and bin size provided
    double bin_edge=min_bin;
    while (bin_edge<max_bin) {
      _bins.push_back(bin_edge);
      bin_edge += bin_size;
    }
    if (abs(_bins.back()-max_bin)>1e-10) _bins.push_back(max_bin);
    
    // initialize values with zeros. Note that _values is one element longer than _bins, since we keep a bin both for underflow (leftmost) and for overflow (rightmost)               
    _values.resize(_bins.size()+1);
    std::fill(_values.begin(),_values.end(),0.0);

    _outfile.open(filename.c_str());
    assert(_outfile.is_open() && "cannot open the specified file.");
  }

  void fill(double xvalue, double yvalue);
  void add(vector<double> additional_histogram);
  void add(vector<double> additional_histogram, double weight);
  void multiply(double multiplier);
  void write();

  int getbin(double value);

  void set_to_zero();

  vector<double> _bins;
  vector<double> _values;
  ofstream _outfile;
};

#endif // __HISTOGRAMS_HH__  
