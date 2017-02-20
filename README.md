# Introduction

The purpose of this project is to provide a platform for automated discovery of [quantum circuits](https://en.wikipedia.org/wiki/Quantum_circuit) by means of a [genetic evolutionary algorithm](https://en.wikipedia.org/wiki/Genetic_algorithm). Defining a simulated population to consist of initially randomly designed quantum circuits and an environment making those individual circuits providing a result closer to a designated goal more likely to survive and procreate, the selection pressure along with the random alterations (mutations, cross-overs) gradually eliminates unsuitable circuits and retains, improves and combines elements of those circuits which perform well. Ideally, this randomly driven process produces (over the course of many generations) solutions very close to perfect.

Currently two real-world applications are demonstrated in this repository:

1. Quantum search circuit discovery, aiming to evolve towards [Grover's algorithm](https://arxiv.org/abs/quant-ph/9605043),

2. Quantum factoring circuit discovery, aiming to produce the core of [Shor's algorithm](https://arxiv.org/abs/quant-ph/9508027).

However, the purpose of this document is to explain how these two goals (with known results) only illustrate the functionality of the more general platform, and provide instruction for implementing evolutions with custom goals.

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

The project itself requires minimal dependencies besides a C++11-compliant compiler and headers (G++ 4.8.1 and higher, Clang 3.3 and higher) and a POSIX-compliant system (Unix, Linux, Mac). It has only been thoroughly tested on Linux + G++ so any feedback from different configurations is most welcome. Windows platform is currently unsupported for the lack of need; in case of interest please file a request in the [issue tracker](https://github.com/vasekp/quantum-ga/issues).

However, the simulation vitally depends on a quantum simulation backend. Currently two backends are supported,

* QIClib (project page on [GitHub](https://titaschanda.github.io/QIClib/)), based on the [Armadillo](http://arma.sourceforge.net/) linear algebra library,

* Quantum++ ([arXiv](https://arxiv.org/abs/1412.4704), [GitHub](https://github.com/vsoftco/qpp)), based on [Eigen3](http://eigen.tuxfamily.org/).

Please refer to the relevant project pages for instructions how to ensure their respective dependencies are properly installed and set up on your system. The default configuration uses QIClib and further relies on [OpenBLAS](http://www.openblas.net/) implementation of the Armadillo routines to reach maximum speed. On a RPM-based system, the dependencies can be installed through the packages `armadillo-devel` and `openblas-devel`. In order to compile with the Quantum++ backend, please refer to the [Makefile](https://github.com/vasekp/quantum-ga/blob/master/Makefile). Note that Fourier transform is faulty and unsupported in recent releases of Eigen3 (3.2.9, 3.2.10) and temporarily disabled in the backend's source code.

# Project structure

The quantum circuit evolution framework is built atop of a concurrently developed [genetic evolution framework](https://github.com/vasekp/genetic), using the high-level genetic evolution interface to wrap around a quantum circuit simulator to evaluate the fitness of individual candidates, which in turn is backed by one of the independent quantum simulation libraries (see System requirements above). This modular design reduces the main program file [quantum.cpp](https://github.com/vasekp/quantum-ga/blob/master/quantum.cpp) virtually just to the core program logic: main program loop and user interface. The numerical settings controlling the genetic evolution (population size, selection pressure, etc.) can also be adjusted in this file.

The specific problem to be solved by the genetic evolution is defined separately in a header file in [include/QGA_Problem/](https://github.com/vasekp/quantum-ga/tree/master/include/QGA_Problem). This needs to define the `Candidate` class which exposes the fitness function. (A candidate base class containing the most functionality is defined in the common base).



# Running

After successful installation, run `./search` or `./fourier` without further parameters. The genetic evolution will commence immediately, starting with a random population. The evaluation of one generation usually takes several hundredths of seconds so the output flows quickly and only summative information is displayed for an overview of the progress and is colour-coded for clarity. On a more concrete example:

![aa](http://i.imgur.com/PV3Bj5q.png)

* **A** denotes the current generation number,
* **B** shows the population size (after culling redundant copies of equal candidates),
* **C** shows the properties of the best-so-far candidate (nondominated and with minimal error),
* **D** shows the size of the nondominated front (the internal population, or archive),
* **E** displays the newest addition to the front,
* **F** is a text-based visualisation of the candidate summarized in **C**.

Due to the nature of the problem and to the nature of the search, there is no quick and guaranteed way of finding the perfect circuit for a given task. A great part of the configuration space needs to be explored before exploiting the discovered features to approach an optimal solution. For this reason a strategy of taking many small steps in a large number of generations has been chosen over a small number of generations employing very elaborate genetic operations. Thus the search may run for a few thousand generations if a perfect solution is required (typically between 1000 and 2000 for the two benchmark problems on 3 qubits, taking 1 to 2 minutes of run time on a modern 4-core processor).

There is no hard-coded termination condition. Surpassing a given error bound and checking for stalled evolution have been considered and rejected. Instead, the user is given a liberty of interrupting and examining the evolution at any point, and to a very limited extent, direct intervention into the evolution is allowed as well. Also see below for termination of the program.

## Interruptions

To pause or stop the evolution and examine the results, hit `Ctrl+C`. The following menu appears:

```
Computation stopped. Choose action:
a: abort,
c: continue,
d: diagnose / list current results,
e: evaluate a candidate in full,
f: filter the front on fitness,
i: inject a candidate,
l: list 20 random candidates,
p: pretty-print a candidate as a circuit,
r: restart,
t: format a candidate as a LuaLaTeX Q-circuit,
q: quit after this generation.
```

User input is then expected in the form of a single lower-case letter followed by `Enter`. If, at this point, you want to exit the program, use **q** (natural stop) or **a** (forced exit) or `Ctrl+C` for a second time (ditto).

The main go-to option to examine the nondominated front is **d**. This lists all the nondominated candidate circuits sorted by their error from highest to lowest. (A high-error candidate can still be a member of the front when it's nondominated in other fitness aspects, e.g., number of gates. An empty circuit usually appears at the top of the list.) Some summative information follows, like the total number of evaluated candidates and time taken.

The list provided by **d** illustrates the format of textual encoding of candidate solutions that is accepted where a candidate is input explicitly by the user (**e**, **i**, **p**, or **t**). The details depend on the problem but the syntax follows a general pattern:

```
SWAP23 Y1(0.6156π) Y2(-0.0021π) P14(0.2028π) Y4(-0.3333π) SWAP14
```

Listing the gates from left to right, i.e., in the order they are applied on the initial state, first comes the name of each, followed by the 1-based qubit indices it acts on. This can be a single digit optionally followed by control qubit indices enclosed in brackets `[`...`]`, or, if all the affected qubits are treated equally, like in a control-Z gate or a swap gate, just a set of the qubit indices. (This assumes no more than 9 qubit lines will be needed, otherwise the current system needs a redesign.) Finally, if a gate has one or more continuous angle parameters, these appear in round parentheses as a multiple of π. Although this is always formatted with a fixed precision and the symbol for π, neither is required when parsing user input. However, no extra space is permitted.

A final choice that deserves some attention in this manual is **f**. Like **d**, this allows to examine the current front, but with the added benefit of filtering on an upper bound of some fitness aspects. For example, in the search problem, one may be interested only in those solutions which use 2 or less oracle calls and don't surpass an error of `0.5`. One could surely list the whole front and ignore solutions which don't qualify but the **f** choice simplifies this. Given that the latter value comes first and the former third within the fitness vector, a specification of such filter would be

```
0.5 _ 2
```

(Parts of fitness which come after the last one we're interested in can be left out.)

# Extending

# Customizing
