PROGS = quantum
HEADERS = include/*.hpp

# Targets can be specified as:
# make qic - use QIClib (default)
# make qpp - use Quantum++
# make qic.d - use QIClib, turn on debug
# make qpp.d - use Quantum++ with debug
#
# Note that all ultimately lead to "quantum" so Make may decide it's not worth 
# recompiling the target. If the platform or debug flag is the only change 
# made, use -B to force a remake.

TARGETS := qpp qic fourier search
TARGETS += $(addsuffix .d,$(TARGETS))

default: search

qpp qic qpp.d qic.d: $(PROGS) $(HEADERS)

CXXFLAGS += -std=c++11 -march=native
CXXFLAGS += -pedantic -Wall -Wextra #-Weffc++
CXXFLAGS += -fno-diagnostics-show-caret -fopenmp
CXXFLAGS += -Iframework/include

$(filter-out %.d,$(TARGETS)): CXXFLAGS += -O3
$(filter %.d,$(TARGETS)): CXXFLAGS += -O1 -g
$(filter qic%,$(TARGETS)): CXXFLAGS += -Iqiclib/include -lopenblas -DUSE_QICLIB
$(filter fourier%,$(TARGETS)): CXXFLAGS += -Iqiclib/include -lopenblas -DUSE_QICLIB -DFOURIER
$(filter search%,$(TARGETS)): CXXFLAGS += -Iqiclib/include -lopenblas -DUSE_QICLIB -DSEARCH
$(filter qpp%,$(TARGETS)): CXXFLAGS += -isystem /usr/include/eigen3 -Iquantum++/include -fpermissive -DUSE_QPP
qpp: CXXFLAGS += -DNODEBUG -DEIGEN_NO_DEBUG
search fourier qic: CXXFLAGS += -DQICLIB_NO_DEBUG

$(TARGETS): $(PROGS)

clean:
	rm $(PROGS)

.PHONY: $(TARGETS) clean

$(PROGS): %: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@
