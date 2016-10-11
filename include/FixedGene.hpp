namespace QGA {

struct gate_struct {
  Backend::Gate op;
  std::string name;
  int inv;
  int sq;
};

std::vector<gate_struct> gates {
  { Backend::I, "I", 0, 0 },
  { Backend::H, "H", 0, -1 },
/*{ Backend::X, "X", 0, -2 },
  { Backend::Y, "Y", 0, -3 },
  { Backend::Z, "Z", 0, -4 },*/
  { Backend::T, "T", +1, 0/*+2*/ },
  { Backend::Ti, "Ti", -1, 0/*+2*/ },
/*{ Backend::S, "S", +1, -3 },
  { Backend::Si, "Si", -1, -4 }*/
};


template<class GeneBase>
class FixedGene : public GeneBase {

  size_t op;
  unsigned tgt;
  unsigned hw;
  Backend::Controls ixs;

  using SP = std::shared_ptr<GeneBase>;

public:

  static SP getNew() {
    /* Distributions: cheap and safer in MT environment this way */
    // distribution of possible gates (except of identity)
    std::uniform_int_distribution<size_t> dOp{1, gates.size() - 1};
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    // distribution of controls
    unsigned tgt_ = dTgt(gen::rng);
    QGA::controls_distribution dCtrl{Config::nBit, Config::pControl, tgt_};
    return std::make_shared<FixedGene>(
        dOp(gen::rng), tgt_, dCtrl(gen::rng));
  }

  Backend::State applyTo(const Backend::State& psi) const override {
    return psi.apply_ctrl(gates[op].op, ixs, tgt);
  }

  unsigned complexity() const override {
    return hw * hw;
  }

  SP invert(const SP& orig) override {
    int dIx = gates[op].inv;
    if(dIx)
      return std::make_shared<FixedGene>(op + dIx, tgt, hw, ixs);
    else
      return orig;
  }

  SP invite(const SP& g) const override {
    return g.get()->visit(g, *this);
  }

  SP visit(const SP& self, const FixedGene& g) override {
    if(op == 0) {
      // Identity * G = G
      return std::make_shared<FixedGene>(g);
    } else if(g.op == 0) {
      // G * Identity = G
      return self;
    } else if(g.op == op && g.tgt == tgt && g.ixs == ixs && gates[op].sq != 0) {
      // G * G = square(G) if also among our operations
      return std::make_shared<FixedGene>(op + gates[op].sq, tgt, hw, ixs);
    } else
      return self;
  }

  std::ostream& write(std::ostream& os) const override {
    os << gates[op].name << tgt + 1;
    if(ixs.size()) {
      os << '[';
      for(auto& ctrl : ixs.as_vector())
        os << ctrl + 1;
      os << ']';
    }
    return os;
  }

  NOINLINE FixedGene(size_t op_, unsigned tgt_, std::vector<bool> ctrl):
      op(op_), tgt(tgt_), ixs(ctrl) {
    hw = ixs.size();
  }

  FixedGene(size_t op_, unsigned tgt_, unsigned hw_,
      const Backend::Controls& ixs_):
    op(op_), tgt(tgt_), hw(hw_), ixs(ixs_) { }

}; // class Gene

} // namespace QGA
