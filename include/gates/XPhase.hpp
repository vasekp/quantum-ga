namespace QGA {


template<class GateBase>
class X : public GateBase {

  unsigned tgt;
  double angle;
  Backend::Gate mat;

  using SP = std::shared_ptr<GateBase>;

public:

  static SP getNew() {
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    // distribution of angle
    std::uniform_real_distribution<> dAng{-0.5*Const::pi, 0.5*Const::pi};
    return std::make_shared<X>(dTgt(gen::rng), dAng(gen::rng));
  }

  Backend::State applyTo(const Backend::State& psi) const override {
    return psi.apply_ctrl(mat, {}, tgt);
  }

  bool isTrivial() override {
    return angle == 0;
  }

  unsigned complexity() const override {
    return 0;
  }

  void invert(SP& self) override {
    self = std::make_shared<X>(tgt, -angle);
  }

  void mutate(SP& self) override {
    std::normal_distribution<> dAng{0.0, 0.1};
    self = std::make_shared<X>(tgt, angle + dAng(gen::rng));
  }

  void simplify(SP& self) override {
    self = std::make_shared<X>(tgt, Tools::rationalize_angle(angle));
  }

  bool invite(SP& first, SP& second) const override {
    return first->merge(first, second, *this);
  }

  bool merge(SP& first, SP& /*second*/, const X& g) override {
    if(g.tgt == tgt) {
      first = std::make_shared<X>(tgt, angle + g.angle);
      return true;
    } else
      return false;
  }

  std::ostream& write(std::ostream& os) const override {
    return os << 'X' << tgt + 1 << '(' << angle / Const::pi << "π)";
  }

  NOINLINE X(unsigned tgt_, double angle_): tgt(tgt_), angle(angle_) {
    mat = Backend::xrot(angle);
  }

}; // class X


template<class GateBase>
class CPhase : public GateBase {

  unsigned tgt;
  double angle;
  Backend::Controls ixs;
  Backend::Gate mat;

  using SP = std::shared_ptr<GateBase>;

public:

  static SP getNew() {
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    // distribution of controls
    unsigned tgt_ = dTgt(gen::rng);
    Tools::controls_distribution dCtrl{Config::nBit, Config::pControl, tgt_};
    // distribution of angle
    std::uniform_real_distribution<> dAng{-0.5*Const::pi, 0.5*Const::pi};
    // Convert P2[13] to P1[23]: mathematically identical and more easily
    // mergeable
    std::vector<bool> ctrl = dCtrl(gen::rng);
    ctrl[tgt_] = true;
    {
      unsigned i;
      for(i = 0; i < Config::nBit; i++)
        if(ctrl[i])
          break;
      // i guaranteed to be < nBit now (at least one bit is set)
      tgt_ = i;
    }
    ctrl[tgt_] = false;
    return std::make_shared<CPhase>(tgt_, dAng(gen::rng), ctrl);
  }

  Backend::State applyTo(const Backend::State& psi) const override {
    return psi.apply_ctrl(mat, ixs, tgt);
  }

  unsigned complexity() const override {
    return ixs.size() * ixs.size();
  }

  void invert(SP& self) override {
    self = std::make_shared<CPhase>(tgt, -angle, ixs);
  }

  void mutate(SP& self) override {
    std::normal_distribution<> dAng{0.0, 0.1};
    self = std::make_shared<CPhase>(tgt, angle + dAng(gen::rng), ixs);
  }

  void simplify(SP& self) override {
    self = std::make_shared<CPhase>(tgt, Tools::rationalize_angle(angle), ixs);
  }

  bool invite(SP& first, SP& second) const override {
    return first->merge(first, second, *this);
  }

  bool merge(SP& first, SP& /*second*/, const CPhase& g) override {
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

  NOINLINE CPhase(unsigned tgt_, double angle_, std::vector<bool> ctrl):
      tgt(tgt_), angle(angle_), ixs(ctrl) {
    mat = Backend::zrot(angle);
  }

  NOINLINE CPhase(unsigned tgt_, double angle_, const Backend::Controls& ixs_):
      tgt(tgt_), angle(angle_), ixs(ixs_) {
    mat = Backend::zrot(angle);
  }

}; // class CPhase

} // namespace QGA
