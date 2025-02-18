# pythia 8 path
PYTHIA8_INCLUDE=/mnt/users/brewerj/pythia8311/include
PYTHIA8_LIB=/mnt/users/brewerj/pythia8311/lib

# fastjet compilation
FASTJET_INCLUDE=/mnt/users/brewerj/fastjet-install/include
FASTJET_LIB=/mnt/users/brewerj/fastjet-install/lib

CXX = g++
CXXFLAGS = -Iinclude -I${PYTHIA8_INCLUDE} -I${FASTJET_INCLUDE} -w -O2 -std=c++17 -pedantic -W -Wall -Wextra -Wshadow -fPIC -pthread

LDFLAGS = -L${PYTHIA8_LIB} -Wl,-rpath,${PYTHIA8_LIB} -L${FASTJET_LIB} -Wl,-rpath,${FASTJET_LIB}
LDLIBS = -lpythia8 -ldl -lfastjet -lRecursiveTools -lFlavorCone -lfastjettools -lgsl

SRCS = src/global_event_analysis.cc src/ccbar_analysis.cc src/EEC.cc src/splitting.cc src/histograms.cc src/medium_mod.cc src/complex_Ei.cc src/jet_energy_loss.cc
OBJS = $(SRCS:src/%.cc=build/%.o)

# Define the two executables and their respective main files
TARGETS = build/compute_substructure build/compute_EEC build/compute_spectra build/simple_energy_loss
MAIN_SUB = src/compute_substructure.cc
MAIN_EEC = src/compute_EEC.cc
MAIN_SPEC = src/compute_spectra.cc
MAIN_SIMPLE = src/simple_energy_loss.cc

# Default target - build both executables
all: $(TARGETS)

build/compute_substructure: $(OBJS) build/compute_substructure.o
	$(CXX) $(OBJS) build/compute_substructure.o -o build/compute_substructure $(LDFLAGS) $(LDLIBS)

build/compute_EEC: $(OBJS) build/compute_EEC.o
	$(CXX) $(OBJS) build/compute_EEC.o -o build/compute_EEC $(LDFLAGS) $(LDLIBS)

build/compute_spectra: $(OBJS) build/compute_spectra.o
	$(CXX) $(OBJS) build/compute_spectra.o -o build/compute_spectra $(LDFLAGS) $(LDLIBS)

build/simple_energy_loss: $(OBJS) build/simple_energy_loss.o
	$(CXX) $(OBJS) build/simple_energy_loss.o -o build/simple_energy_loss $(LDFLAGS) $(LDLIBS)

build/compute_substructure.o: $(MAIN_SUB)
	$(CXX) $(CXXFLAGS) -c $(MAIN_SUB) -o build/compute_substructure.o

build/compute_EEC.o: $(MAIN_EEC)
	$(CXX) $(CXXFLAGS) -c $(MAIN_EEC) -o build/compute_EEC.o

build/compute_spectra.o: $(MAIN_SPEC)
	$(CXX) $(CXXFLAGS) -c $(MAIN_SPEC) -o build/compute_spectra.o

build/simple_energy_loss.o: $(MAIN_SIMPLE)
	$(CXX) $(CXXFLAGS) -c $(MAIN_SIMPLE) -o build/simple_energy_loss.o

# Rule to compile the shared modules
build/%.o: src/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f build/*.o $(TARGETS)
