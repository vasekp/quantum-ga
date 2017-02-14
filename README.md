# Introduction

The purpose of this project is to provide a platform for automated discovery of [quantum circuits]
(https://en.wikipedia.org/wiki/Quantum_circuit) by means of a [genetic evolutionary algorithm]
(https://en.wikipedia.org/wiki/Genetic_algorithm). Defining a simulated population to consist of initially randomly designed
quantum circuits and an environment making those individual circuits providing a result closer to a designated goal more likely
to survive and procreate, the selection pressure along with the random alterations (mutations, cross-overs) gradually eliminates
unsuitable circuits and retains, improves and combines elements of those circuits which perform well. Ideally, this randomly
driven process produces (over the course of many generations) solutions very close to perfect.

Currently two real-world applications are demonstrated in this repository:

1. Quantum search circuit discovery, aiming to evolve towards [Grover's algorithm](https://arxiv.org/abs/quant-ph/9605043),

2. Quantum factoring circuit discovery, aiming to produce the core of [Shor's algorithm](https://arxiv.org/abs/quant-ph/9508027).

However, the purpose of this document is to explain how these two goals (with known results) only illustrate the functionality
of the more general platform, and provide instruction for implementing evolutions with custom goals.

# Installation

This project is distributed in source code form. To download, use

```
git clone --recurse-submodules https://github.com/vasekp/quantum-ga
```

To compile, use

```
make all
```

To run (Grover Search), use

```
./search
```

## System requirements

The project itself requires minimal dependencies besides a C++11-compliant compiler and headers (G++ 4.8.1 and higher, Clang
3.3 and higher) and a POSIX-compliant system (Unix, Linux, Mac). Windows platform is currently unsupported for the lack of
need; in case of interest please file a request in the [issue tracker](https://github.com/vasekp/quantum-ga/issues).

However, the simulation vitally depends on a quantum simulation backend. Currently two backends are supported,

* QIClib (project page on [GitHub](https://titaschanda.github.io/QIClib/)), based on the
[Armadillo](http://arma.sourceforge.net/) linear algebra library,

* Quantum++ ([arXiv](https://arxiv.org/abs/1412.4704), [GitHub](https://github.com/vsoftco/qpp)), based on
[Eigen3](http://eigen.tuxfamily.org/).

Please refer to the relevant project pages for instructions how to ensure their respective dependencies are properly installed
and set up on your system. The default configuration uses QIClib and further relies on [OpenBLAS](http://www.openblas.net/)
implementation of the Armadillo routines to reach maximum speed. On a RPM-based system, the dependencies can be installed through
the packages `armadillo-devel` and `openblas-devel`. In order to compile with the Quantum++ backend, please refer to the
[Makefile](https://github.com/vasekp/quantum-ga/blob/master/Makefile).

# Project structure

The quantum circuit evolution framework is built atop of a concurrently developed [genetic evolution framework]
(https://github.com/vasekp/genetic), using the high-level genetic evolution interface to wrap around a quantum circuit simulator
to evaluate the fitness of individual candidates, which in turn is backed by one of the independent quantum simulation libraries
(see System requirements above).
This modular design reduces the main program file [quantum.cpp](https://github.com/vasekp/quantum-ga/blob/master/quantum.cpp)
virtually just to the core program logic: main program loop and user interface. The numerical settings controlling the genetic
evolution (population size, selection pressure, etc.) can also be adjusted in this file.

The specific problem to be solved by the genetic evolution is defined separately in a header file in [include/QGA_Problem/]
(https://github.com/vasekp/quantum-ga/tree/master/include/QGA_Problem). This needs to define the `Candidate` class which
exposes the fitness function. (A candidate base class containing the most functionality is defined in the common base).
