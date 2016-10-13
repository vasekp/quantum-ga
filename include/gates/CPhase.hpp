#ifndef GATE_CPHASE_HPP
#define GATE_CPHASE_HPP

namespace QGA {

namespace Gates {

using Tools::Controls;

template<class GateBase>
class CPhase : public GateBase {

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
    Tools::controls_distribution<Controls::ANY>
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
    return std::make_shared<CPhase>(tgt, dAng(gen::rng), ctrl);
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

  void invert(Pointer& self) const override {
    self = std::make_shared<CPhase>(tgt, -angle, ixs);
  }

  void mutate(Pointer& self) const override {
    std::normal_distribution<> dAng{0.0, 0.1};
    self = std::make_shared<CPhase>(tgt, angle + dAng(gen::rng), ixs);
  }

  void simplify(Pointer& self) const override {
    self = std::make_shared<CPhase>(tgt, Tools::rationalize_angle(angle), ixs);
  }

  bool invite(Pointer& first, Pointer& second) const override {
    return first->merge(first, second, *this);
  }

  bool merge(Pointer& first, Pointer&, const CPhase& g) const override {
    if(g.tgt == tgt && g.ixs == ixs) {
      first = std::make_shared<CPhase>(tgt, angle + g.angle, ixs);
      return true;
    } else
      return false;
  }

  std::ostream& write(std::ostream& os) const override {
    os << "P[" << tgt + 1;
    for(auto ctrl : ixs.as_vector())
      os << ctrl + 1;
    os << "](" << angle / Const::pi << "π)";
    return os;
  }

  CPhase(unsigned tgt_, double angle_, std::vector<bool> ctrl):
      tgt(tgt_), angle(angle_), ixs(ctrl), mat(Backend::phase(angle)) { }

  CPhase(unsigned tgt_, double angle_, const Backend::Controls& ixs_):
      tgt(tgt_), angle(angle_), ixs(ixs_), mat(Backend::phase(angle)) { }

}; // class CPhase

} // namespace Gates

} // namespace QGA

#endif // !defined GATE_CPHASE_HPP
