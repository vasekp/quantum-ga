namespace QGA {

namespace Gates {

struct CPhase {

template<class GateBase>
class Inner : public GateBase {

  unsigned tgt;
  double angle;
  Backend::Controls ixs;
  Backend::Gate mat;

  using typename GateBase::Pointer;
  using Ctx = typename GateBase::Context;

public:

  // construct a random gate
  Inner():
    tgt(std::uniform_int_distribution<unsigned>{0, Config::nBit - 1}
        (gen::rng)),
    angle(angle_distribution<>{}(gen::rng)),
    ixs(), mat(func::phase(angle))
  {
    // distribution of controls
    controls_distribution<Controls::ANY>
      dCtrl{Config::nBit, tgt, Config::pControl};
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
  Inner(unsigned tgt_, double angle_, const Backend::Controls& ixs_):
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

  Pointer invert(const Pointer&) const override {
    return std::make_shared<Inner>(tgt, -angle, ixs);
  }

  Pointer mutate(const Pointer&) const override {
    angle_distribution<true> dAng{};
    return std::make_shared<Inner>(tgt, angle + dAng(gen::rng), ixs);
  }

  Pointer simplify(const Pointer&) const override {
    return std::make_shared<Inner>(tgt, rationalize_angle(angle), ixs);
  }

  void hit(typename GateBase::Counter& c) const {
    c.hit(this);
  }

  Pointer invite(const Pointer& first) const override {
    return first->merge(*this);
  }

  Pointer merge(const Inner& g) const override {
    if(g.tgt == tgt && g.ixs == ixs)
      return std::make_shared<Inner>(tgt, angle + g.angle, ixs);
    else
      return {};
  }

  std::ostream& write(std::ostream& os) const override {
    os << "P" << tgt + 1;
    for(auto ctrl : ixs.as_vector())
      os << ctrl + 1;
    os << "(" << angle / Const::pi << "π)";
    return os;
  }

  static Pointer read(const std::string& s) {
    std::string reS{};
    std::regex re{"P(\\d+)\\((-?[0-9.]+)(π)?\\)"};
    std::smatch m{};
    if(!std::regex_match(s, m, re))
      return {};
    std::vector<bool> ctrl(Config::nBit, false);
    unsigned tgt = (unsigned)(~0);
    for(const char& c : m[1].str()) {
      size_t pos = c - '1';
      if(pos >= 0 && pos < Config::nBit && pos != tgt) {
        if(tgt == (unsigned)(~0))
          tgt = pos;
        else
          ctrl[pos] = true;
      }
    }
    double angle = std::stod(m[2].str()) * Const::pi;
    return std::make_shared<Inner>(tgt, angle, Backend::Controls{ctrl});
  }

}; // class CPhase::Inner

template<class GateBase>
using Template = Inner<GateBase>;

}; // struct CPhase

} // namespace Gates

} // namespace QGA
