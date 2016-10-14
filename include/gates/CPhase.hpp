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

public:

  static Pointer getNew() {
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    // distribution of controls
    unsigned tgt = dTgt(gen::rng);
    controls_distribution<Controls::ANY>
      dCtrl{Config::nBit, tgt, Config::pControl};
    // distribution of angle
    std::uniform_real_distribution<> dAng{-0.5*Const::pi, 0.5*Const::pi};
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
    return std::make_shared<Inner>(tgt, dAng(gen::rng), ctrl);
  }

  Backend::State applyTo(const Backend::State& psi) const override {
    return psi.apply_ctrl(mat, ixs, tgt);
  }

  bool isTrivial() const override {
    return angle == 0;
  }

  unsigned complexity() const override {
    return ixs.size() * ixs.size();
  }

  Pointer invert(const Pointer&) const override {
    return std::make_shared<Inner>(tgt, -angle, ixs);
  }

  Pointer mutate(const Pointer&) const override {
    std::normal_distribution<> dAng{0.0, 0.1};
    return std::make_shared<Inner>(tgt, angle + dAng(gen::rng), ixs);
  }

  Pointer simplify(const Pointer&) const override {
    return std::make_shared<Inner>(tgt, rationalize_angle(angle), ixs);
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
    os << "P[" << tgt + 1;
    for(auto ctrl : ixs.as_vector())
      os << ctrl + 1;
    os << "](" << angle / Const::pi << "π)";
    return os;
  }

  Inner(unsigned tgt_, double angle_, const Backend::Controls& ixs_):
      tgt(tgt_), angle(angle_), ixs(ixs_), mat(Backend::phase(angle)) { }

}; // class CPhase::Inner

template<class GateBase>
using Template = Inner<GateBase>;

}; // struct CPhase

} // namespace Gates

} // namespace QGA
