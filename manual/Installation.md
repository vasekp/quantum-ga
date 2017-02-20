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

Please refer to the relevant project pages for instructions how to ensure their respective dependencies are properly installed and set up on your system. The default configuration uses QIClib and further relies on [OpenBLAS](http://www.openblas.net/) implementation of the Armadillo routines to reach maximum speed. On a RPM-based system, the dependencies can be installed through the packages `armadillo-devel` and `openblas-devel`. In order to compile with the Quantum++ backend, please consult the [Makefile](https://github.com/vasekp/quantum-ga/blob/master/Makefile). Note that Fourier transform is faulty and unsupported in recent releases of Eigen3 (3.2.9, 3.2.10) and temporarily disabled in the backend's source code.

- - -

Back to [the README](https://github.com/vasekp/quantum-ga/blob/readme/README.md)
