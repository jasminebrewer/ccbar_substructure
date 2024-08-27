# Substructure and energy correlators of doubly-charm tagged jets

This repository provides the code necessary to study the properties of the $g \to c\bar{c}$ splitting function in vacuum, and to compute its modification in the quark-gluon plasma. This repository contains code for two example analyses in this framework
1. \underline{c-cbar jet substructure}: identify jets containing a c-cbar pair, reconstruct the features of the $g \to c\bar{c}$ splitting according to Pythia, and use jet reclustering to provide a data-driven reconstruction of the features of the splitting. Compute the medium modification of the splitting.
2. \underline{energy correlators of c-cbar jets}: identify jets containing a c-cbar pair, reconstruct the features of the $g \to c\bar{c}$ splitting according to Pythia, use tagged charm quarks (or hadrons) to compute two- and three-point energy correlators.

## Getting Started:
This repository requires that Pythia8 and fastjet are installed in order to run. The user should edit `PYTHIA8_LIB`, `PYTHIA8_INCLUDE`, `FASTJET_LIB`, and `FASTJET_INCLUDE` in the Makefile to point to their local installation.
From the main directory, then create a build directory and build the executables
```
mkdir build
make
```
This builds executables both to computer c-cbar substructure (`build/compute_substructure`) and to compute the 2- and 3-point energy correlators of c-cbar jets (`build/compute_EEC`)

Provided with the repository are an example parameter file, `example_params.dat`, and an example shell script `run.sh`. The meaning of each parameter in the parameter file is discussed further below. `run.sh` takes as command line arguments 'substructure' or 'EEC' (depending on which analysis should be run), and the parameter file for the run, and saves the output of the run in the directory runs/(substructure or EEC)/run_*, which * is a 6-digit code distinguishing each run.

### Explanation:
The medium modification is implemented as described in (...)

### Parameters:
`maxneventsperjob`: number of Pythia events to generate

`tune`: Pythia tune to use

`beamidA`, `beamidB`: particle ids of the incoming beams (2212 for protons)

`eCM`: center of mass energy of the collision, in GeV

`switchoffbhadrondecays`: turns off the decays of B hadrons (only applicable when hadronization is on).

`pTHatMin`: minimum transverse momentum (pT) of the event generation. For physical results, should be less than the minimum jet pT.

`ISR`: flag for initial state radiation (allowed values: on and off)

`MPI`: flag for multiple-parton interactions (allowed values: on and off)

`FSR`: flag for final-state radiation (allowed values: on and off)

`jetR`: radius to use for jet finding

`trackPtMin`: minimum allowed transverse momentum for tracks (particles)

`trackEtaCut`: maximum allowed absolute value pseudorapidity for tracks (particles)

`HFPtMin`: minimum allowed transverse momentum for heavy-flavor particles (note that "heavy flavor" is defined by eventSelection, so if eventSelection is "qq" then this is the minimum tranvserse momentum for the light quarks specified by the event selection

`jetPtMin`: minimum transverse momentum of jets

`jetEtaMin`, `jetEtaMax`: minimum and maximum values for the pseudorapidity of jets

`jetAlgorithm`: jet algorithm to use for the jet finding. Allowed values are antikt, kt, and CA

`jetReclAlgorithm`: jet algorithm to use for the jet reclustering step (after the jet finding). Allowed values are antikt, kt, and CA

`process`: specified the allowed processes in the Pythia hard event event generation. Allowed values are "all" (HardQCD:all=on) or "gg" (gluon initial states only)

`eventSelection`: flag to change the event selection criteria for parton-level events. Allowed values are "qq" (light quarks, specifically down), "cc" (charm quarks), and "bb" (bottom quarks)

`hadronization`: flag to specify if the event generation should be parton-level or hadron-level. Allowed values are on and off.

`qhatL`: product of the medium parameter qhat and the length of the medium L. Used for the medium modification.

`L`: length of the medium, in fermi

Leaving out parameters for the Pythia generation will revert to Pythia defaults.
