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
class XYZGene : public GeneBase {

  size_t op;
  double angle;
  double gphase;
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
    std::uniform_int_distribution<unsigned> dCtrl{};
    // distribution of angle
    std::uniform_real_distribution<> dAng{-0.5*Const::pi, 0.5*Const::pi};
    return std::make_shared<XYZGene>(
        dOp(gen::rng), dAng(gen::rng), dAng(gen::rng),
        dTgt(gen::rng), dCtrl(gen::rng));
  }

  Backend::State applyTo(const Backend::State& psi) const override {
    return psi.apply_ctrl(mat, ixs, tgt);
  }

  unsigned complexity() const override {
    return hw * hw;
  }

  /*double phase() const {
    return gphase;
  }*/

  SP invert(const SP&) override {
    return std::make_shared<XYZGene>(op, -angle, -gphase, tgt, ixs, hw);
  }

  SP mutate(const SP&) override {
    std::normal_distribution<> dAng{0.0, 0.1};
    std::bernoulli_distribution dWhich{};
    bool modAngle = dWhich(gen::rng);
    return std::make_shared<XYZGene>(
        op,
        modAngle ? angle + dAng(gen::rng) : angle,
        modAngle ? gphase : gphase + dAng(gen::rng),
        tgt, ixs, hw);
  }

  SP simplify(const SP&) override {
    return std::make_shared<XYZGene>(op,
        rationalize(std::fmod(angle / Const::pi, 2.0)) * Const::pi,
        rationalize(std::fmod(gphase / Const::pi, 2.0)) * Const::pi,
        tgt, ixs, hw);
  }

  SP invite(const SP& g) const override {
    return g.get()->visit(g, *this);
  }

  SP visit(const SP& self, const XYZGene& g) override {
    if(angle == 0) {
      // op1 = (phase*)identity
      return std::make_shared<XYZGene>(
          g.op, g.angle, g.gphase + gphase, g.tgt, g.ixs, g.hw);
    } else if(g.angle == 0) {
      // op2 = (phase*)identity
      return std::make_shared<XYZGene>(
          op, angle, gphase + g.gphase, tgt, ixs, hw);
    } else if(g.op == op && g.tgt == tgt && g.ixs == ixs) {
      return std::make_shared<XYZGene>(
          op, angle + g.angle, gphase + g.gphase, tgt, ixs, hw);
    } else
      return self;
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

  NOINLINE XYZGene(size_t op_, double angle_, double phase_,
    unsigned tgt_, unsigned control_enc):
      op(op_), angle(angle_), gphase(phase_), tgt(tgt_), hw(0), ixs() {
    if(gates[op].ctrl) {
      std::vector<bool> bits{GeneBase::ctrlBitString(control_enc, tgt)};
      ixs = Backend::Controls{bits};
      hw = ixs.size();
    }
    computeMat();
  }

  NOINLINE XYZGene(size_t op_, double angle_, double phase_,
    unsigned tgt_, const Backend::Controls& ixs_, unsigned hw_):
      op(op_), angle(angle_), gphase(phase_), tgt(tgt_), hw(hw_), ixs(ixs_) {
    computeMat();
  }

private:

  double rationalize(double x) {
    double a = std::abs(x);
    constexpr unsigned N = 8;
    double coeffs[N];
    unsigned t;
    for(t = 0; t < N; t++) {
      coeffs[t] = std::floor(a);
      if(coeffs[t] > 100) {
        coeffs[t++] = 100;
        break;
      }
      a = 1/(a - coeffs[t]);
    }
    std::discrete_distribution<unsigned> dStop(&coeffs[1], &coeffs[t]);
    unsigned cut = dStop(gen::rng) + 1;
    if(cut == t)
      return x;
    a = coeffs[--cut];
    while(cut > 0)
      a = coeffs[--cut] + 1/a;
    return x < 0 ? -a : a;
  }

  void computeMat() {
    if(op >= gates.size())
      throw std::logic_error("gate must be between 0 and 2");
    mat = std::exp(gphase * Const::i) * gates[op].fn(angle);
  }

}; // class XYZGene

} // namespace QGA
