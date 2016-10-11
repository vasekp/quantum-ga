namespace QGA {

struct gate_struct {
  Backend::Gate(*fn)(double);
  std::string name;
  bool ctrl;
};

std::vector<gate_struct> gates {
  {Backend::xrot, "X", true},
  {Backend::yrot, "Y", true},
  {Backend::zrot, "Z", true}
};


template<class GeneBase>
class XYZ : public GeneBase {

  size_t op;
  double angle;
  unsigned tgt;
  unsigned hw;
  Backend::Controls ixs;
  Backend::Gate mat;

  using SP = std::shared_ptr<GeneBase>;

public:

  static SP getNew() {
    // distribution of possible gates
    std::uniform_int_distribution<size_t> dOp{0, gates.size() - 1};
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    // distribution of controls
    unsigned tgt_ = dTgt(gen::rng);
    QGA::controls_distribution dCtrl{Config::nBit, Config::pControl, tgt_};
    // distribution of angle
    std::uniform_real_distribution<> dAng{-0.5*Const::pi, 0.5*Const::pi};
    return std::make_shared<XYZ>(
        dOp(gen::rng), dAng(gen::rng), tgt_, dCtrl(gen::rng));
  }

  Backend::State applyTo(const Backend::State& psi) const override {
    return psi.apply_ctrl(mat, ixs, tgt);
  }

  unsigned complexity() const override {
    return hw * hw;
  }

  void invert(SP& self) override {
    self = std::make_shared<XYZ>(op, -angle, tgt, ixs, hw);
  }

  void mutate(SP& self) override {
    std::normal_distribution<> dAng{0.0, 0.1};
    self = std::make_shared<XYZ>(op, angle + dAng(gen::rng), tgt, ixs, hw);
  }

  void simplify(SP& self) override {
    self = std::make_shared<XYZ>(op,
        GeneBase::rationalize(std::fmod(angle / Const::pi, 2.0)) * Const::pi,
        tgt, ixs, hw);
  }

  bool invite(SP& first, SP& second) const override {
    return first->visit(first, second, *this);
  }

  bool visit(SP& first, SP& second, const XYZ& g) override {
    if(angle == 0) {
      // op1 = identity: replace by second and consume
      first = second;
      return true;
    } else if(g.angle == 0) {
      // op2 = identity: consume
      return true;
    } else if(g.op == op && g.tgt == tgt && g.ixs == ixs) {
      first = std::make_shared<XYZ>(op, angle + g.angle, tgt, ixs, hw);
      return true;
    } else
      return false;
  }

  std::ostream& write(std::ostream& os) const override {
    os << gates[op].name << tgt + 1;
    if(ixs.size()) {
      os << '[';
      for(auto ctrl : ixs.as_vector())
        os << ctrl + 1;
      os << ']';
    }
    os << '(' << angle / Const::pi << "π)";
    return os;
  }

// Should be private, but constructors are needed my std::make_shared().

  NOINLINE XYZ(size_t op_, double angle_, unsigned tgt_,
      std::vector<bool> ctrl):
      op(op_), angle(angle_), tgt(tgt_), hw(0), ixs(ctrl) {
    if(gates[op].ctrl)
      hw = ixs.size();
    else
      ixs.clear();
    mat = gates[op].fn(angle);
  }

  NOINLINE XYZ(size_t op_, double angle_, unsigned tgt_,
      const Backend::Controls& ixs_, unsigned hw_):
      op(op_), angle(angle_), tgt(tgt_), hw(hw_), ixs(ixs_) {
    mat = gates[op].fn(angle);
  }

}; // class XYZ

} // namespace QGA
