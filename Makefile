PROGS = quantum
HEADERS = include/*.hpp

all: CXXFLAGS += -O3
debug: CXXFLAGS += -O1 -g

all debug: $(PROGS)

CXXFLAGS += -std=c++11 -march=native
CXXFLAGS += -pedantic -Wall -Wextra #-Weffc++
CXXFLAGS += -fno-diagnostics-show-caret -fopenmp
CXXFLAGS += -Iframework/include
#CXXFLAGS += -isystem /usr/include/eigen3 -Iquantum++/include -DNODEBUG -DEIGEN_NO_DEBUG 
CXXFLAGS += -Iqiclib/include -lopenblas

clean:
	rm $(PROGS)

.PHONY: all debug clean

$(PROGS): %: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@
