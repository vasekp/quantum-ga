# Extending

The system has been designed with extensibility in mind. The primary goal is that the two example problems, Fourier and Search, will serve as templates on how to build other tasks for the genetic search. The latter in addition uses a nonstandard custom gate (the oracle) and as such can be consulted on the practical details of the implementation thereof.

## Step-by-step specification of a problem

The task for which the framework should seek a quantum circuit is specified via a `Candidate` class. The [main program](https://github.com/vasekp/quantum-ga/blob/master/quantum.cpp) references this class but its implementation is let unspecified, expected to be declared in an included header. However, a base class template `QGA::CandidateBase` is provided by the framework headers providing most of the functionality so that really only problem-specific components need to be provided.

The two example problems are both found under [include/QGA_Problem](https://github.com/vasekp/quantum-ga/tree/master/include/QGA_Problem), although this is no fixed requirement. Let us walk through the implementation of the Fourier problem to see a particular program in action.

```c++
#ifndef QGA_PROBLEM_HPP
#define QGA_PROBLEM_HPP

namespace {
```
These lines should appear at the top of each problem file.

```c++
  using Gene = QGA::Gene<QGA::Gates::Y, QGA::Gates::CPhase, QGA::Gates::SWAP>;
```
The `QGA::Gene` template represents objects participating in the genotype, i.e., individual circuit gates. This is where the [gate types](#gate-types) that are allowed in the evolution are specified. The way the `QGA::Gene` class template is designed it needs to be aware of its descendants (for the purposes of static polymorphism).

```c++
  class Candidate : public QGA::CandidateBase<Candidate, Gene, double, double> {
```
The `Candidate` class is the centre of the definition of the Fourier problem. It must derive from `QGA::Candidate`, specifying as template parameters
* itself (for [CRTP](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern)),
* the gene class (containing in it the gate types), as declared above,
* a flat list of the fitness element types, in this case average error and maximum error, both `double`. Note that this will automatically be augmented by gate counts by type (i.e. count of Y-rotations, count of controlled phases, and swap count, in this order).

```c++
  using Base = QGA::CandidateBase<Candidate, Gene, double, double>;

public:
  Base::FitnessMain fitness_main() const {
    if(genotype().size() > 1000)
      return {INFINITY, INFINITY};
    using cxd = std::complex<double>;
    cxd overlapTotal{0};
    double errorMax = 0;
    unsigned dim = 1 << Config::nBit;
    State psi{};
    for(unsigned i = 0; i < dim; i++) {
      psi.reset(i);
      State out = State::fourier(psi);
      cxd overlap = State::overlap(out, sim(psi));
      overlapTotal += overlap;
      double error = 1.0 - std::abs(overlap);
      if(error > errorMax)
        errorMax = error;
    }
    double errorAvg = std::max(1.0 - std::abs(overlapTotal / cxd(dim)), 0.0);
    return {
      this->trimError(errorAvg),
      this->trimError(errorMax)
    };
  }
```
This is the fitness evaluation function. Using the base class's `FitnessMain` return type, which encapsulates two values of type `double` (as requested in the specialization of `QGA::CandidateBase`), it calculates the "average error" and "maximum error" and returns them as a tuple. Note that the naming is arbitrary in the sense that the former is not an average of a sample of which the latter would be maximum: in this case there is some additional logic designed to accept solutions which represent the correct unitary transform up to a global complex phase. In this logic, an overlap is calculated as a complex number for the image of each basis state to the desired output state. An error is defined as the complement of the absolute of an overlap value into 1. The average error is then simply the error of an averaged overlap as a complex number. (This is somewhat deceptive because an average error can be larger than the maximum error, e.g., if each of the transformed vectors is a multiple of the desired output but the relative phases are wrong.) Also note that the function immediately rejects some candidates based on an upper limit of their genotype size. This is a basic measure to prevent [bloat](http://cswww.essex.ac.uk/staff/poli/gp-field-guide/113Bloat.html) and to reject solutions which achieve some imaginary improvement by accumulating numerical errors (as sometimes indeed happens).

The functions `State::overlap` or `State::reset`, as well as the constructor, are examples of the abstraction provided by the backend layer. These map (usually rather trivially) to functions of the underlying quantum simulation or linear algebra packages.

The fitness function calls a private class member function `sim` which is omitted for brevity. It can be found in [Fourier.hpp](https://github.com/vasekp/quantum-ga/blob/master/include/QGA_Problem/Fourier.hpp). Finally,

```c++
  std::ostream& print_full(std::ostream& os) const {
    unsigned dim = 1 << Config::nBit;
    State psi{};
    os << '\n';
    for(unsigned i = 0; i < dim; i++) {
      psi.reset(i);
      State out = sim(psi);
      for(unsigned j = 0; j < dim; j++)
        os << std::abs(out[j])*std::sqrt(dim) << "/√" << dim << "∠"
          << std::showpos << std::arg(out[j]) / QGA::Const::pi << "π "
          << std::noshowpos;
      os << '\n';
    }
    return os;
  }
```
is a function called when a full listing of a candidate's results is required. This happens at the exit of the program when a perfect or almost perfect solution is found, or when the program is [interrupted](https://github.com/vasekp/quantum-ga/blob/master/manual/Running.md#interruptions) with the **e** command (evaluate a candidate in full). This is just a callback of the standard `ostream` output operator and shares its signature. The candidate is quite free in implementation. This version outputs the transform of each basis vector using polar representation of each amplitude probability.

## Gate types

The following is a list of the currently implemented quantum gate types intended for generic usage. They are all defined in the [include/QGA_bits/gates/](https://github.com/vasekp/quantum-ga/tree/master/include/QGA_bits/gates/) directory.
* `QGA::XYZ`, `QGA::X`, `QGA::Y`, `QGA::Z`: parametric Pauli rotations,
* `QGA::Fixed`: fixed 1-qubit gate set (the identity, Hadamard, X, Y, Z Pauli gates, T and S phase rotations and their inverses),
* `QGA::SU2`: generic SU(2) matrix (beware: expands the configuration space drastically and tends to converge slowly),
* `QGA::CNOT`: Controlled-NOT 2-qubit gate,
* `QGA::SWAP`: 2-qubit swap,
* `QGA::CPhase`: symmetric *k*-qubit controlled phase gate.

Most of these gates come in a generic setting that can be further adjusted. This is done via an inner class alias mechanism. For example, to allow control qubits for phase rotations, replace `QGA::XYZ` by
```c++
QGA::XYZ::WithControls<QGA::Controls::ANY>
```
The `WithControls` specifier can be used with all the above gate types except `QGA::SWAP` and takes the following constants of the `QGA::Controls` enum:
* `NONE`: no control bits allowed (this is the default for `XYZ` and derivates, `Fixed`, `SU2`),
* `ONE`: exactly one control bit required (this is the default for `CNOT`),
* `ANY`: any combination of control bits, including zero control bits, is allowed (this is the default for `CPhase`),
* `LEAST1`: any combination of control bits but at least one must be present.

Similar to `WithControls` there is `WithGates` which is the mechanisms through which the generic `QGA::XYZ` is internally restricted to `QGA::X` and others. This can be used, for example, when introducing a new parametric rotation gate or when designing a custom gate set for `QGA::Fixed`. For the details, see the source of [Fixed.hpp](https://github.com/vasekp/quantum-ga/blob/master/include/QGA_bits/gates/Fixed.hpp).

The specifiers can also be chained and even repeated, with the last occurrence of each overriding the previous instances of the same.

When these gate types are not sufficient, it is possible to extend the list by defining custom ones. For details please see the corresponding section of the [Customizing](https://github.com/vasekp/quantum-ga/blob/master/manual/Customizing.md#defining-custom-gate-types) guide.

- - -

Back to [the README](https://github.com/vasekp/quantum-ga/blob/master/README.md)
