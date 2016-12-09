namespace QGA {

namespace Gates {

namespace internal {

template<Controls cc>
class CNOT {

template<class GateBase>
class CNOTTemp : public GateBase {

  using typename GateBase::Pointer;
  using Ctx = typename GateBase::Context;

public:

  // construct a random gate
  CNOTTemp():
    tgt(std::uniform_int_distribution<unsigned>{0, Config::nBit - 1}
        (gen::rng)),
    ixs(controls_distribution<cc>{Config::nBit, tgt, Config::pControl}
        (gen::rng)),
    odd(true)
  { }

  // construct using parameters
  CNOTTemp(unsigned tgt_, const Backend::Controls& ixs_, bool odd_ = true):
      tgt(tgt_), ixs(ixs_), odd(odd_) { }

  CNOTTemp(bool): tgt(), ixs(), odd(false) { }

  Backend::State applyTo(const Backend::State& psi, const Ctx*) const override {
    return odd ? psi.apply_ctrl(Backend::X, ixs, tgt) : psi;
  }

  bool isTrivial() const override {
    // CNOT^(2k) = CNOT^0 = identity
    return !odd;
  }

  unsigned controls() const override {
    return ixs.size();
  }

  Pointer mutate(const Pointer&) const override {
    return std::make_shared<CNOTTemp>();
  }

  Pointer swapQubits(const Pointer& self, unsigned s1, unsigned s2)
    const override
  {
    return odd
      ? std::make_shared<CNOTTemp>(
          tgt == s1 ? s2 : tgt == s2 ? s1 : tgt,
          Backend::Controls::swapQubits(ixs, s1, s2))
      : self;
  }

  void hit(typename GateBase::Counter& c) const {
    c.hit(this);
  }

  const CNOTTemp* cast(const CNOTTemp*) const override {
    return this;
  }

  bool sameType(const GateBase& other) const override {
    const CNOTTemp* c = other.cast(this);
    return c != nullptr && c->tgt == tgt && c->ixs == ixs;
  }

  Pointer merge(const GateBase& other) const override {
    if(!sameType(other))
      return {};
    const CNOTTemp* c = other.cast(this);
    return std::make_shared<CNOTTemp>(tgt, ixs, odd ^ c->odd);
  }

  std::ostream& write(std::ostream& os) const override {
    if(!odd)
      return os << "[Id]";
    os << "NOT" << tgt + 1;
    if(ixs.size()) {
      os << '[';
      for(auto ctrl : ixs.as_vector())
        os << ctrl + 1;
      os << ']';
    }
    return os;
  }

  void printOn(QGA::CircuitPrinter& p) const override {
    if(odd)
      p.addControlledGate("X", tgt, ixs.as_vector());
  }

  static Pointer read(const std::string& s) {
    regex::regex re{"(\\[Id\\])|NOT(\\d)(\\[(\\d+)\\])?"};
    regex::matches ms{};
    if(!re.match(s, ms))
      return {};
    if(ms.matched(1))
      return std::make_shared<CNOTTemp>(false);
    unsigned tgt = ms.match(2)[0] - '1';
    if(tgt >= Config::nBit)
      return {};
    std::vector<bool> ctrl(Config::nBit, false);
    if(ms.matched(3))
      for(const char& c : ms.match(4)) {
        size_t pos = c - '1';
        if(pos >= 0 && pos < Config::nBit && pos != tgt)
          ctrl[c - '1'] = true;
      }
    return std::make_shared<CNOTTemp>(tgt, Backend::Controls{ctrl});
  }

private:

  unsigned tgt;
  Backend::Controls ixs;
  bool odd;  // parity of the power

}; // class CNOT<Controls>::CNOTTemp<GateBase>

template<class GateBase>
using Template = CNOTTemp<GateBase>;

template<Controls cc_>
using WithControls = CNOT<cc_>;

}; // struct CNOT<Controls>

} // namespace internal

using CNOT = internal::CNOT<Controls::ONE>;

} // namespace Gates

} // namespace QGA
