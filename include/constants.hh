#ifndef __CONSTANTS_HH__
#define __CONSTANTS_HH__

namespace constants
{
  // identifying particle ids used by pythia
  const int CHARM = 4;
  const int ANTICHARM = -4;
  const int BOTTOM = 5;
  const int ANTIBOTTOM = -5;
  const int DOWN = 1;
  const int ANTIDOWN = -1;
  const int GLUON = 21;
  const int D0 = 421;
  const int D0BAR = -421;
  const int B0 = 511;
  const int B0BAR = -511;
  // particles masses
  const double CHARM_MASS = 1.27;
  const double BOTTOM_MASS = 4.18;
  // status code for particles in pythia
  const int INCOMING_HARD = -21; // incoming from hard process
  const int INCOMING_SUB = -31; // incoming from subprocess
  const int INCOMING_ISR = -41; // incoming from initial state radiation
  const int OUTGOING_HARD = -23; // outgoing from hard process
  const int OUTGOING_SUB = -33; // outgoing from subprocess
  const int OUTGOING_ISR = -43; // outgoing from initial state radiation
  const int SHOWER = -51; // produced by the shower
  // unit conversions
  const double invGeVtofm = 0.1973; // conversion for length units in GeV^{-1} to fermi
  // QCD constants
  const double CA = 3.;
  const double CF = 4./3.;
  const double Nc = 3.;
}

#endif // __CONSTANTS_HH__