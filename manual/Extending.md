# Extending

The system has been designed with extensibility in mind. The primary goal is that the two example problems, Fourier and Search, will serve as templates on how to build other tasks for the genetic search. The latter in addition uses a nonstandard custom gate (the oracle) and as such can be consulted on the practical details of the implementation thereof.

## Step-by-step specification of a problem

The task for which the framework should seek a quantum circuit is specified via a `Candidate` class. The [main program](https://github.com/vasekp/quantum-ga/blob/master/quantum.cpp) references this class but its implementation is let unspecified, expected to be declared in an included header. However, a base class template `QGA::CandidateBase` is provided by the framework headers providing most of the functionality so that really only problem-specific components need to be provided.

The two example problems are both found under [include/QGA_Problem](https://github.com/vasekp/quantum-ga/tree/master/include/QGA_Problem), although this is no fixed requirement. Let us walk through the implementation of the Fourier problem to see a particular program in action.

```
#ifndef QGA_PROBLEM_HPP
#define QGA_PROBLEM_HPP

namespace {
```
These lines should appear at the top of each problem file.

```
using Gene = QGA::Gene<QGA::Gates::Y, QGA::Gates::CPhase, QGA::Gates::SWAP>;
```
The `QGA::Gene` template represents objects participating in the genotype, i.e., individual circuit gates. This is where the [gate types](#gate-types) that are allowed in the evolution are specified. The way the `QGA::Gene` class template is designed it needs to be aware of its descendants (for the purposes of static polymorphism).

```
class Candidate : public QGA::CandidateBase<Candidate, Gene, double, double> {
```
The `Candidate` class is the centre of the definition of the Fourier problem. It must derive from `QGA::Candidate`, specifying as template parameters
* itself (for [CRTP](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern)),
* the gene class (containing in it the gate types), as declared above,
* a flat list of the fitness element types, in this case average error and maximum error, both `double`. Note that this will automatically be augmented by gate counts by type (i.e. count of Y-rotations, count of controlled phases, and swap count, in this order).

```
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

```
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
is a function called when a full listing of a candidate's results is required. This happens at the exit of the program when a perfect or almost perfect solution is found, or when the program is [interrupted](https://github.com/vasekp/quantum-ga/blob/readme/manual/Running.md#interruptions) with the **e** command (evaluate a candidate in full). This is just a callback of the standard `ostream` output operator and shares its signature. The candidate is quite free in implementation. This version outputs the transform of each basis vector using polar representation of each amplitude probability.

## Gate types

The following is a list of the currently implemented quantum gate types intended for generic usage. They are all defined in the [include/QGA_bits/gates/](https://github.com/vasekp/quantum-ga/tree/master/include/QGA_bits/gates/) directory.
* `QGA::XYZ`, `QGA::X`, `QGA::Y`, `QGA::Z`: parametric Pauli rotations,
* `QGA::Fixed`: fixed 1-qubit gate set (the identity, Hadamard, X, Y, Z Pauli gates, T and S phase rotations and their inverses),
* `QGA::SU2`: generic SU(2) matrix (beware: expands the configuration space drastically and tends to converge slowly),
* `QGA::CNOT`: Controlled-NOT 2-qubit gate,
* `QGA::SWAP`: 2-qubit swap,
* `QGA::CPhase`: symmetric *k*-qubit controlled phase gate.

Most of these gates come in a generic setting that can, however, be further adjusted. This is done via an inner class alias mechanism. For example, to allow control qubits for phase rotations, replace `QGA::XYZ` by
```
QGA::XYZ::WithControls<QGA::Controls::ANY>
```
The `WithControls` specifier can be used with all the above gate types except `QGA::SWAP` and takes the following constants of the `QGA::Controls` enum:
* `NONE`: no control bits allowed (this is the default for `XYZ` and derivates, `Fixed`, `SU2`),
* `ONE`: exactly one control bit required (this is the default for `CNOT`),
* `ANY`: any combination of control bits, including zero control bits, is allowed (this is the default for `CPhase`),
* `LEAST1`: any combination of control bits but at least one must be present.

Similar to `WithControls` there is `WithGates` which is the mechanisms through which the generic `QGA::XYZ` is internally restricted to `QGA::X` and others. This can be used, for example, when introducing a new parametric rotation gate or when designing a custom gate set for `QGA::Fixed`. For the details, see the source of [Fixed.hpp](https://github.com/vasekp/quantum-ga/blob/master/include/QGA_bits/gates/Fixed.hpp).

The specifiers can also be chained and even repeated, with the last occurrence of each overriding the previous instances of the same.

## Defining custom gate types

Sometimes a problem requires a very specific gate type which is not covered by the base set but would not be beneficial enough to other simulations. This is the case of the oracle gate of the search problem: it takes an external parameter which needs to be adjusted between runs on an otherwise identical input and circuit, and different oracle gates within a single circuit need to be linked to the same external parameter. It also acts on all qubits at once (NB: an oracle affecting the relative phase of the marked node is implemented, which does not rely on an ancilla qubit) and thus behaves quite differently from the other gates. On the other hand, it does not have any internal configuration like rotation angles (by its logic, it is a black box element): only its presence matters. Such gate can be defined directly in the problem file and then fed among others to the `QGA::CandidateBase` template just as described above. Let us study this on the case of [Search.hpp](https://github.com/vasekp/quantum-ga/blob/master/include/QGA_Problem/Search.hpp).

```
struct Oracle {

template<class GateBase>
class OracleTemp : public GateBase {
```
Gates are somewhat trickier to implement due to the design internal to the QGA framework: each "gate type" is in fact just a wrapper for a template which receives a base class as its parameter. This, again, is a necessity due to the adopted static polymorphism pattern used to massively speed up the genetic operations in the presence of principally different gate types. The template is expected to be called `Template` but a more descriptive name has been chosen here, with an alias near the end of the `Oracle` class.

```
  using typename GateBase::Pointer;
  using typename GateBase::Context;
```
These two class aliases provide some insight into the internals of the framework's design. `Pointer` is an indirect pointer to an instance of a gate; in fact a `std::shared_ptr` to `GateBase`. It is the mechanism through which member functions of a gate are passed instances of other gates where needed (e.g., for comparison). `Context` wraps any circuit-wide information that needs to be passed to its constituing gates. This is usually just `void` and ignored but it's very important in the case of the search problem where the `Context` is a structure containing `unsigned mark`, the marked node.

```
  bool odd;  // parity of the power
```
This line is interesting. We will see references to the value of `odd` several times below. This is a switch on the action of the oracle. If the power is even (`false`), signifying an even number of oracle gates in sequence, these undo each other's action and the result is identity. An oracle gate is by default generated with an odd power. However this mechanism becomes used when neighbouring oracle gates are merged into one.

```
public:
  OracleTemp(bool odd_ = true): odd(odd_) { }
```
The sole constructor of the `OracleTemp` class template.

```
  State applyTo(const State& psi, const Context* pMark) const override {
    unsigned mark = pMark->mark;
    State ret{psi};
    if(odd)
      ret[mark] = -ret[mark];
    return ret;
  }
```
This is the most important function of the `Oracle` gate: here the gate gets applied to a state vector, returning by value. We extract the `mark` value from the passed `Context` and use this to apply a relative phase to the corresponding base vector component of `psi`. (This happens when the power is odd, otherwise `psi` is returned unchanged.) The bracket access to the state vector's element is another abstraction unified by the backend.

```
  bool isTrivial() const override {
    // oracle^(2k) = oracle^0 = identity
    return !odd;
  }
```
This function, needed by `merge` below, returns whether a given instance represents a trivial (identity) gate. Such gates are usually artifacts of special cases of previous merges and are short-lived; they are typically not found in the output. The oracle becomes an identity when the power is even.

```
  Pointer swapQubits(const Pointer& self, unsigned, unsigned) const override {
    return self;
  }
```
The function `swapQubits` is called when two input qubits are to be exchanged on a gate. Its usual task is to reroute control or target qubits. Given that an oracle gate uses all qubit lines and treats them equally, its instance can be returned unchanged. There is one more important observation happening here. Note that the return reference happens through the `Pointer` class, aliased in the top of the class. A self-pointer is passed as an argument, too; returning `this` would break `std::shared_ptr`'s reference counting mechanism.

The two following functions need to appear verbatim in each gate:
```
  void hit(typename GateBase::Counter& c) const {
    c.hit(this);
  }

  const OracleTemp* cast(const OracleTemp*) const override {
    return this;
  }
```
The first is responsible for gate counting (for the purposes of fitness evaluation) and can not appear in the base class because it resolves on the derived class type of `this`. The latter is an implementation of the static polymorphism pattern.

The following two functions are used in gate merging mechanisms:
```
  bool sameType(const GateBase& other) const override {
    const OracleTemp* c = other.cast(this);
    return c != nullptr && c->odd == odd;
  }

  Pointer merge(const GateBase& other) const override {
    if(!sameType(other))
      return {};
    // oracle * oracle = oracle^2 → true ^ true = false
    const OracleTemp* c = other.cast(this);
    return std::make_shared<OracleTemp>(odd ^ c->odd);
  }
```
Their names are self-documenting but their implementation deserves some attention. Namely, the static up-cast provided by the `GateBase::cast` function. Here we pass `this` to request the target type. If the `other` gate is an `OracleTemp`, the overload defined above is called, returning its pointer. For all other derived classes, the default behaviour is to return a `nullptr`. So, in order to check whether `this` and `other` are of the same type, we check the result against a `nullptr`. Moreover, all trivial gates are considered a separate category (identity), so we (conditioned on the non-`nullptr` check) compare these data members, too.

`merge` attempts to merge a gate with another, which is to be applied *after* it. The return type is a `Pointer` but neither the self-pointer nor the `Pointer` to the `other` gate are provided. This is because the cases where trivially one gate or the other can be returned through their `Pointer` are sorted out before control is given to the `merge` overload, thus its function is in most cases to create an entirely new gate or to report a failure. The latter is done by calling the default constructor of `Pointer`, as in the third line. A new gate is created via a call to `std::make_shared`. In this case, we call the only constructor `OracleTemp::OracleTemp(bool = true)`.

What is left are functions relating to textual output and input:
```
  std::ostream& write(std::ostream& os) const override {
    return os << (odd ? "Oracle" : "[Id]");
  }

  void printOn(QGA::CircuitPrinter& p) const override {
    if(odd)
      p.addBarrierGate("U_f");
  }

  static Pointer read(const std::string& s) {
    regex::regex re{"\\[Id\\]|(Oracle)"};
    regex::matches ms{};
    if(!re.match(s, ms))
      return {};
    return std::make_shared<OracleTemp>(ms.matched(1));
  }
```
The `write` and `read` functions extend the mechanism described in [Interruptions](https://github.com/vasekp/quantum-ga/blob/readme/manual/Running.md#interruptions). The output of `write` is required to be parsed correctly by `read`, and should never collide with the description of any other gate type, the only permitted exception is `[Id]` for trivial gates. For the description of the circuit-formatting functions of `QGA::CircuitPrinter` see the interface defined in [its source code](https://github.com/vasekp/quantum-ga/blob/master/include/CircuitPrinter.hpp).

Besides `merge` there are other logical functions that a gate type can define, with their default implementations returning a no-op (where possible) or failing safely. See `invert` and `simplify` in [GateBase.hpp](https://github.com/vasekp/quantum-ga/blob/master/include/QGA_bits/GateBase.hpp) for details.

This concludes the definition of `Oracle::OracleTemp` template, but the outer class still needs to define the last alias
```
template<class GateBase>
using Template = OracleTemp<GateBase>;
```
Note that this would not be needed would `OracleTemp` be directly called `Template`, but the former name sheds more logic onto what's happening in the code. Also, other gates use this space to define the class aliases like `WithControls` etc.

- - -

Back to [the README](https://github.com/vasekp/quantum-ga/blob/readme/README.md)
