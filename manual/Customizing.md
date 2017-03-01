# Customizing

## Controlling the evolution parameters

The genetic evolution (population size, archive size, selection pressure, etc.) is controlled by a set of parameters declared in [QGA_commons.hpp](https://github.com/vasekp/quantum-ga/blob/master/include/QGA_commons.hpp) and defined near the beginning of [the main program file](https://github.com/vasekp/quantum-ga/blob/master/quantum.cpp). Please refer to the comments in the latter file with regard to customizing their values.

All of these currently need to be specified as compile-time constants. This has been a convenience when designing the framework, but for all save a small number there is no reason preventing their input via command line, or changing them at runtime, in a future version. This should be very straightforward to change if needed in a particular application. The only exception where the code does not except a runtime change of value is `Config::nBit`, the number of qubits in the system. Lifting this restriction will be the immediate focus of near-future development.

## Defining custom gate types

Sometimes a problem requires a very specific gate type which is not covered by the base set but would not be beneficial enough to other simulations. This is the case of the oracle gate of the search problem: it takes an external parameter which needs to be adjusted between runs on an otherwise identical input and circuit, and different oracle gates within a single circuit need to be linked to the same external parameter. It also acts on all qubits at once (NB: an oracle affecting the relative phase of the marked node is implemented, which does not rely on an ancilla qubit) and thus behaves quite differently from the other gates. On the other hand, it does not have any internal configuration like rotation angles (by its logic, it is a black box element): only its presence matters. Such gate can be defined directly in the problem file and then fed among others to the `QGA::CandidateBase` template just as described above. Let us study this on the case of [Search.hpp](https://github.com/vasekp/quantum-ga/blob/master/include/QGA_Problem/Search.hpp).

```c++
struct Oracle {

template<class GateBase>
class OracleTemp : public GateBase {
```
Gates are somewhat trickier to implement due to the design internal to the QGA framework: each "gate type" is in fact just a wrapper for a template which receives a base class as its parameter. This, again, is a necessity due to the adopted static polymorphism pattern used to massively speed up the genetic operations in the presence of principally different gate types. The template is expected to be called `Template` but a more descriptive name has been chosen here, with an alias near the end of the `Oracle` class.

```c++
  using typename GateBase::Pointer;
  using typename GateBase::Context;
```
These two class aliases provide some insight into the internals of the framework's design. `Pointer` is an indirect pointer to an instance of a gate; in fact a `std::shared_ptr` to `GateBase`. It is the mechanism through which member functions of a gate are passed instances of other gates where needed (e.g., for comparison). `Context` wraps any circuit-wide information that needs to be passed to its constituing gates. This is usually just `void` and ignored but it's very important in the case of the search problem where the `Context` is a structure containing `unsigned mark`, the marked node.

```c++
  bool odd;  // parity of the power
```
This line is interesting. We will see references to the value of `odd` several times below. This is a switch on the action of the oracle. If the power is even (`false`), signifying an even number of oracle gates in sequence, these undo each other's action and the result is identity. An oracle gate is by default generated with an odd power. However this mechanism becomes used when neighbouring oracle gates are merged into one.

```c++
public:
  OracleTemp(bool odd_ = true): odd(odd_) { }
```
The sole constructor of the `OracleTemp` class template. In general there can be more, but a default constructor (no parameters) must be available. If the gate has some internal parameters the default constructor should initialize them randomly.

```c++
  State applyTo(const State& psi, const Context* pMark) const override {
    unsigned mark = pMark->mark;
    State ret{psi};
    if(odd)
      ret[mark] = -ret[mark];
    return ret;
  }
```
This is the most important function of the `Oracle` gate: here the gate gets applied to a state vector, returning by value. We extract the `mark` value from the passed `Context` and use this to apply a relative phase to the corresponding base vector component of `psi`. (This happens when the power is odd, otherwise `psi` is returned unchanged.) The bracket access to the state vector's element is another abstraction unified by the backend.

```c++
  bool isTrivial() const override {
    // oracle^(2k) = oracle^0 = identity
    return !odd;
  }
```
This function, needed by `merge` below, returns whether a given instance represents a trivial (identity) gate. Such gates are usually artifacts of special cases of previous merges and are short-lived; they are typically not found in the output. The oracle becomes an identity when the power is even.

```c++
  Pointer swapQubits(const Pointer& self, unsigned, unsigned) const override {
    return self;
  }
```
The function `swapQubits` is called when two input qubits are to be exchanged on a gate. Its usual task is to reroute control or target qubits. Given that an oracle gate uses all qubit lines and treats them equally, its instance can be returned unchanged. There is one more important observation happening here. Note that the return reference happens through the `Pointer` class, aliased in the top of the class. A self-pointer is passed as an argument, too; returning `this` would break `std::shared_ptr`'s reference counting mechanism.

The two following functions need to appear verbatim in each gate:
```c++
  void hit(typename GateBase::Counter& c) const {
    c.hit(this);
  }

  const OracleTemp* cast(const OracleTemp*) const override {
    return this;
  }
```
The first is responsible for gate counting (for the purposes of fitness evaluation) and can not appear in the base class because it resolves on the derived class type of `this`. The latter is an implementation of the static polymorphism pattern.

The following two functions are used in gate merging mechanisms:
```c++
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
```c++
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
```c++
template<class GateBase>
using Template = OracleTemp<GateBase>;
```
Note that this would not be needed would `OracleTemp` be directly called `Template`, but the former name sheds more logic onto what's happening in the code. Also, other gates use this space to define the class aliases like `WithControls` etc.

## Customizing genetic operators

All the genetic operator logic is provided by [CandidateFactory.hpp](https://github.com/vasekp/quantum-ga/blob/master/include/QGA_bits/CandidateFactory.hpp). This class provides functionality for selecting candidates from the population and mutating and combining them, as well as a static interface for generating the initial population. Some adaptive heuristics for choosing genetic operators with priority based on their prior success have also been implemented but deprecated since, and may be removed in a future revision.

Wherever a random count of gates is needed, a geometric distribution is used. This is achieved by probing a Bernoulli-distributed boolean random value until first failure. This way a scale-free distribution can be obtained while controlling the expected count remains easy. Depending on the context, the latter is specified by constants `Config::expLengthIni` and `Config::expMutationCount` of [QGA_commons.hpp](https://github.com/vasekp/quantum-ga/blob/master/include/QGA_commons.hpp).

### Enabling and disabling genetic operators

Each genetic operator is represented by one private member function of the `CandidateFactory` class. They are chosen at random in `CandidateFactory::getNew()` by means of the ancillary `Selector` class. The list of enabled operators, along with their display names, appears near the very end of [CandidateFactory.hpp](https://github.com/vasekp/quantum-ga/blob/master/include/QGA_bits/CandidateFactory.hpp).

### Programming a custom genetic operator

Any additional genetic operators should be added to the `CandidateFactory` class's private section and follow the logic of the original set. Let us study, for example, `CandidateFactory::mDeleteSlice`:

```c++
  Candidate mDeleteSlice() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();  
    if(sz == 0)
      return parent;
```
We start with drawing a random candidate from the population using a common mechanism provided by the method `get()`. This and its genotype should be stored by reference for speed. Most genetic operators are conditioned by a minumum size of the genotype (here 1), a fast rejection can be achieved by returning `parent`. (This results in a copy being taken.)

```c++
    std::geometric_distribution<size_t> dGeom{1.0 / Config::expMutationCount};
    std::uniform_int_distribution<size_t> dPos{0, sz - 1};
    size_t pos1 = dPos(gen::rng),        
           len = 1 + dGeom(gen::rng),
           pos2 = pos1 + len > sz ? sz : pos1 + len;
```
Deleting a slice means choosing two nearby points within the genotype and removing what is between them, but at least one gene. We prepare the two markers `pos1` and `pos2` accordingly.

```c++
    std::vector<Gene> gtNew{};       
    gtNew.reserve(sz - (pos2 - pos1));  
    gtNew.insert(gtNew.end(), gtOrig.begin(), gtOrig.begin() + pos1);
    gtNew.insert(gtNew.end(), gtOrig.begin() + pos2, gtOrig.end());
    return Candidate{std::move(gtNew)};
  }
```
Finally, a new genotype is prepared as an empty `std::vector` with a pre-allocated size, and the parts before and after the slice are copied to it. Note that copying genes is a quick operation, as they internally are shared pointers. The candidate to be returned by the operator is prepared using this vector, removing its contents in turn.

### Customizing the selection strategy

The model for population has been chosen to be a NSGA-selection based one. This affects the main program file as well as `CandidateFactory::get()`. If intending to replace this by a different selection mechanism please refer to the underlying genetic framework documentation (available separately [here](https://vasekp.github.io/genetic/doc/index.html)).

- - -

Back to [the README](https://github.com/vasekp/quantum-ga/blob/readme/README.md)
