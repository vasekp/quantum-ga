# To compile against Quantum++ (default is QIClib) use:
# make BACKEND=QPP target
#
# To make a debug version:
# make DEBUG=1 target
#
# To force remake a target (with different defines):
# make -B target

SOURCES = quantum.cpp
HEADERS = include/*.hpp include/*/*.hpp

TARGETS := simple fourier search
default: search

CXXFLAGS += -std=c++11 -march=native
CXXFLAGS += -pedantic -Wall -Wextra -Weffc++
CXXFLAGS += -fopenmp
CXXFLAGS += -Iframework/include
CXXFLAGS += -fno-diagnostics-show-caret -fmax-errors=3

ifdef DEBUG
	CXXFLAGS += -O0 -g
else
	CXXFLAGS += -O3 -DNODEBUG -DEIGEN_NO_DEBUG -DQICLIB_NO_DEBUG
endif

ifeq ($(BACKEND), QPP)
	CXXFLAGS += -isystem /usr/include/eigen3 -Iquantum++/include -DUSE_QPP
else
	# QIClib is the default
	CXXFLAGS += -Iqiclib/include -lopenblas -DUSE_QICLIB
endif

fourier: CXXFLAGS += -DFOURIER

search:	CXXFLAGS += -DSEARCH

all: $(TARGETS)

$(TARGETS): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $@

clean:
	-rm $(TARGETS)

.PHONY: all clean

