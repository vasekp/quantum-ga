# To compile against Quantum++ (default is QIClib) use:
# make BACKEND=QPP target
#
# To make a debug version:
# make DEBUG=1 target
#
# To make a benchmark (a perfectly reproducible run with a fixed thread count 
# and fixed seed of random number generators):
# make BENCH=1 target
#
# To force remake a target (with different defines):
# make touch target

SOURCES = quantum.cpp
HEADERS = include/*.hpp include/*/*.hpp include/*/*/*.hpp
LIBS = regex.o

TARGETS := simple fourier search
default: search

CXXFLAGS += -std=c++11 -march=native
CXXFLAGS += -pedantic -Wall -Wextra -Weffc++
CXXFLAGS += -fopenmp
CXXFLAGS += -Iframework/include -Iinclude
CXXFLAGS += -fno-diagnostics-show-caret -fmax-errors=3

ifdef DEBUG
	CXXFLAGS += -O0 -g
else
	CXXFLAGS += -O3 -DNODEBUG -DEIGEN_NO_DEBUG -DQICLIB_NO_DEBUG
endif

ifdef BENCH
	CXXFLAGS += -DBENCH
endif

ifeq ($(BACKEND), QPP)
	CXXFLAGS += -isystem /usr/include/eigen3 -Iquantum++/include -DUSE_QPP
else
	# QIClib is the default
	CXXFLAGS += -Iqiclib/include -lopenblas -DUSE_QICLIB
endif

ifdef NBIT
	CXXFLAGS += -DNBIT=${NBIT}
endif

fourier: CXXFLAGS += -DFOURIER

search:	CXXFLAGS += -DSEARCH

all: $(TARGETS)

$(TARGETS): $(SOURCES) $(HEADERS) $(LIBS)
	$(CXX) $(CXXFLAGS) $(SOURCES) $(LIBS) -o $@

touch:
	touch $(SOURCES)

regex.o: regex.cpp include/regex.hpp
	$(CXX) -O3 -c regex.cpp

clean:
	-rm $(TARGETS) $(LIBS)

.PHONY: all clean touch
