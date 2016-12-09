namespace QGA {

namespace Gates {

struct gate_struct_f {
  const Backend::Gate* op;
  std::string name;
  int inv;
  int sq;
};


namespace internal {

static const std::vector<gate_struct_f> gates_fixed {
  { &Backend::I, "I", 0, 0 },
  { &Backend::H, "H", 0, -1 },
  { &Backend::X, "X", 0, -2 },
  { &Backend::Y, "Y", 0, -3 },
  { &Backend::Z, "Z", 0, -4 },
  { &Backend::T, "T", +1, +2 },
  { &Backend::Ti, "Ti", -1, +2 },
  { &Backend::S, "S", +1, -3 },
  { &Backend::Si, "Si", -1, -4 }
};

template<Controls cc, const std::vector<gate_struct_f>* gates>
struct Fixed {

template<class GateBase>
class FixedTemp : public GateBase {

  using typename GateBase::Pointer;
  using Ctx = typename GateBase::Context;

public:

  // construct a random gate
  FixedTemp():
    op(std::uniform_int_distribution<size_t>{1, gates->size() - 1}(gen::rng)),
    tgt(std::uniform_int_distribution<unsigned>{0, Config::nBit - 1}
        (gen::rng)),
    ixs(controls_distribution<cc>{Config::nBit, tgt, Config::pControl}
        (gen::rng))
  { }

  // construct using parameters
  FixedTemp(size_t op_, unsigned tgt_, const Backend::Controls& ixs_):
    op(op_), tgt(tgt_), ixs(ixs_) { }

  Backend::State applyTo(const Backend::State& psi, const Ctx*) const override {
    return psi.apply_ctrl(*(*gates)[op].op, ixs, tgt);
  }

  bool isTrivial() const override {
    return op == 0;
  }

  unsigned controls() const override {
    return ixs.size();
  }

  Pointer invert(const Pointer& self) const override {
    int dIx = (*gates)[op].inv;
    if(dIx != 0)
      return std::make_shared<FixedTemp>(op + dIx, tgt, ixs);
    else
      return self;
  }

  Pointer mutate(const Pointer&) const override {
    return std::make_shared<FixedTemp>();
  }

  Pointer swapQubits(const Pointer&, unsigned s1, unsigned s2) const override {
    return std::make_shared<FixedTemp>(
        op,
        tgt == s1 ? s2 : tgt == s2 ? s1 : tgt,
        Backend::Controls::swapQubits(ixs, s1, s2));
  }

  void hit(typename GateBase::Counter& c) const {
    c.hit(this);
  }

  const FixedTemp* cast(const FixedTemp*) const override {
    return this;
  }

  bool sameType(const GateBase& other) const override {
    const FixedTemp* c = other.cast(this);
    return c != nullptr && c->tgt == tgt && c->ixs == ixs && c->op == op;
  }

  Pointer merge(const GateBase& other) const override {
    if(!sameType(other))
      return {};
    // G * G = square(G) if also among our operations
    if((*gates)[op].sq != 0)
      return std::make_shared<FixedTemp>(op + (*gates)[op].sq, tgt, ixs);
    else
      return {};
  }

  std::ostream& write(std::ostream& os) const override {
    os << (*gates)[op].name << tgt + 1;
    if(ixs.size()) {
      os << '[';
      for(auto ctrl : ixs.as_vector())
        os << ctrl + 1;
      os << ']';
    }
    return os;
  }

  void printOn(QGA::CircuitPrinter& p) const override {
    p.addControlledGate((*gates)[op].name, tgt, ixs.as_vector());
  }

  static Pointer read(const std::string& s) {
    std::string reS{};
    for(const gate_struct_f& g : *gates)
      reS = reS + "|(" + g.name + ")";
    regex::regex re{"(?:" + reS.substr(1) + ")(\\d)(\\[(\\d+)\\])?"};
    regex::matches ms{};
    if(!re.match(s, ms))
      return {};
    size_t num = gates->size();
    size_t op;
    for(op = 0; op < num; op++)
      if(ms.matched(op + 1))
        break;
    unsigned tgt = ms.match(num + 1)[0] - '1';
    if(tgt >= Config::nBit)
      return {};
    std::vector<bool> ctrl(Config::nBit, false);
    if(ms.matched(num + 2))
      for(const char& c : ms.match(num + 3)) {
        size_t pos = c - '1';
        if(pos >= 0 && pos < Config::nBit && pos != tgt)
          ctrl[pos] = true;
      }
    return std::make_shared<FixedTemp>(op, tgt, Backend::Controls{ctrl});
  }

private:

  size_t op;
  unsigned tgt;
  Backend::Controls ixs;

}; // class Fixed<Controls, Gates>::FixedTemp<GateBase>

template<class GateBase>
using Template = FixedTemp<GateBase>;

template<Controls cc_>
using WithControls = Fixed<cc_, gates>;

template<const std::vector<gate_struct_f>* gates_>
using WithGates = Fixed<cc, gates_>;

}; // struct Fixed<Controls, Gates>

} // namespace internal

using Fixed = internal::Fixed<Controls::NONE, &internal::gates_fixed>;

} // namespace Gates

} // namespace QGA
