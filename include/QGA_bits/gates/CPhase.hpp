namespace QGA {

namespace Gates {

namespace internal {

template<Controls cc>
struct CPhase {

template<class GateBase>
class CPhaseTemp : public GateBase {

  using typename GateBase::Pointer;
  using Ctx = typename GateBase::Context;

public:

  // construct a random gate
  CPhaseTemp():
    tgt(std::uniform_int_distribution<unsigned>{0, Config::nBit - 1}
        (gen::rng)),
    angle(angle_distribution<>{}(gen::rng)),
    ixs(), mat(func::phase(angle))
  {
    // distribution of controls
    controls_distribution<cc> dCtrl{Config::nBit, tgt, Config::pControl};
    // Convert P2[13] to P1[23]: mathematically identical and more easily
    // mergeable
    std::vector<bool> ctrl = dCtrl(gen::rng);
    ctrl[tgt] = true;
    {
      unsigned i;
      for(i = 0; i < Config::nBit; i++)
        if(ctrl[i])
          break;
      // i guaranteed to be < nBit now (at least one bit was set)
      tgt = i;
    }
    ctrl[tgt] = false;
    ixs = Backend::Controls{ctrl};
  }

  // construct using parameters
  CPhaseTemp(unsigned tgt_, double angle_, const Backend::Controls& ixs_):
      tgt(tgt_), angle(angle_), ixs(ixs_), mat(func::phase(angle)) { }

  Backend::State applyTo(const Backend::State& psi, const Ctx*) const override {
    return psi.apply_ctrl(mat, ixs, tgt);
  }

  bool isTrivial() const override {
    return angle == 0;
  }

  unsigned controls() const override {
    return ixs.size();
  }

  Pointer getAnother() const override {
    return std::make_shared<CPhaseTemp>();
  }

  Pointer invert(const Pointer&) const override {
    return std::make_shared<CPhaseTemp>(tgt, -angle, ixs);
  }

  Pointer mutate(const Pointer&) const override {
    angle_distribution<true> dAng{};
    return std::make_shared<CPhaseTemp>(tgt, angle + dAng(gen::rng), ixs);
  }

  Pointer simplify(const Pointer&) const override {
    return std::make_shared<CPhaseTemp>(tgt, rationalize_angle(angle), ixs);
  }

  Pointer swapQubits(const Pointer& self, unsigned s1, unsigned s2)
    const override
  {
    std::vector<bool> ctrl(Config::nBit);
    for(auto v : ixs.as_vector())
      ctrl[v] = true;
    ctrl[tgt] = true;
    if(ctrl[s1] == ctrl[s2]) // swapping has no effect
      return self;
    // if they are different then swapping amounts to negation
    ctrl[s1] = !ctrl[s1];
    ctrl[s2] = !ctrl[s2];
    unsigned tgt_;
    for(tgt_ = 0; tgt_ < Config::nBit; tgt_++)
      if(ctrl[tgt_])
        break;
    ctrl[tgt_] = false;
    Backend::Controls ixs_{ctrl};
    return std::make_shared<CPhaseTemp>(tgt_, angle, ixs_);
  }

  unsigned type() const override {
    return GateBase::Indexer::index(this);
  }

  const CPhaseTemp* cast(const CPhaseTemp*) const override {
    return this;
  }

  bool sameType(const GateBase& other) const override {
    const CPhaseTemp* c = other.cast(this);
    return c != nullptr && c->tgt == tgt && c->ixs == ixs;
  }

  Pointer merge(const GateBase& other) const override {
    if(!sameType(other))
      return {};
    const CPhaseTemp* c = other.cast(this);
    return std::make_shared<CPhaseTemp>(tgt, angle + c->angle, ixs);
  }

  std::ostream& write(std::ostream& os) const override {
    os << "P" << tgt + 1;
    for(auto ctrl : ixs.as_vector())
      os << ctrl + 1;
    os << "(" << angle / Const::pi << "π)";
    return os;
  }

  void printOn(QGA::CircuitPrinter& p) const override {
    p.addControlledGate("Φ", tgt, ixs.as_vector());
  }

  static Pointer read(const std::string& s) {
    std::string reS{};
    regex::regex re{"P(\\d+)\\((-?[0-9.]+)(π)?\\)"};
    regex::matches ms{};
    if(!re.match(s, ms))
      return {};
    std::vector<bool> ctrl(Config::nBit, false);
    unsigned tgt = (unsigned)(~0);
    for(auto c : ms.match(1)) {
      size_t pos = c - '1';
      if(pos >= 0 && pos < Config::nBit && pos != tgt) {
        if(tgt == (unsigned)(~0))
          tgt = pos;
        else
          ctrl[pos] = true;
      }
    }
    double angle = std::stod(ms.match(2)) * Const::pi;
    return std::make_shared<CPhaseTemp>(tgt, angle, Backend::Controls{ctrl});
  }

private:

  unsigned tgt;
  double angle;
  Backend::Controls ixs;
  Backend::Gate mat;

}; // class CPhase<Controls>::CPhaseTemp<GateBase>

template<class GateBase>
using Template = CPhaseTemp<GateBase>;

template<Controls cc_>
using WithControls = CPhase<cc_>;

}; // struct CPhase<Controls>

} // namespace internal

using CPhase = internal::CPhase<Controls::ANY>;

} // namespace Gates

} // namespace QGA
