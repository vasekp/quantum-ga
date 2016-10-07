namespace Wrapper {

struct gate_struct {
  Gate op;
  std::string name;
  int inv;
  int sq;
};

std::vector<gate_struct> gates {
  { internal::I, "I", 0, 0 },
  { internal::H, "H", 0, -1 },
/*{ internal::X, "X", 0, -2 },
  { internal::Y, "Y", 0, -3 },
  { internal::Z, "Z", 0, -4 },*/
  { internal::T, "T", +1, 0/*+2*/ },
  { internal::Ti, "Ti", -1, 0/*+2*/ },
/*{ internal::S, "S", +1, -3 },
  { internal::Si, "Si", -1, -4 }*/
};


template<class GeneBase>
class FixedGene : public GeneBase {

  size_t op;
  unsigned tgt;
  unsigned hw;
  arma::uvec ixs;

  using SP = std::shared_ptr<GeneBase>;

public:

  static SP getNew() {
    /* Distributions: cheap and safer in MT environment this way */
    // distribution of possible gates
    std::uniform_int_distribution<size_t> dOp{1, gates.size() - 1};
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{1, Config::nBit};
    // distribution of controls
    std::uniform_int_distribution<unsigned> dCtrl{};
    return std::make_shared<FixedGene>(
        dOp(gen::rng), dTgt(gen::rng), dCtrl(gen::rng));
  }

  State apply(const State& psi) const override {
    return qic::apply_ctrl(
        psi,
        gates[op].op,
        ixs,
        {tgt});
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
    } else if(g.op == op
        && g.tgt == tgt
        && g.ixs.size() == ixs.size()
        && arma::all(g.ixs == ixs)
        && gates[op].sq != 0) {
      // G * G = square(G) if also among our operations
      return std::make_shared<FixedGene>(op + gates[op].sq, tgt, hw, ixs);
    } else
      return self;
  }

  std::ostream& write(std::ostream& os) const override {
    os << gates[op].name << tgt;
    if(ixs.size()) {
      os << '[';
      for(auto& ctrl : ixs)
        os << ctrl;
      os << ']';
    }
    return os;
  }

  NOINLINE FixedGene(size_t op_, unsigned tgt_, unsigned control_enc):
      op(op_), tgt(tgt_), hw(0) {
    std::vector<bool> bits{GeneBase::ctrlBitString(control_enc, tgt - 1)};
    std::vector<arma::uword> ixv;
    ixv.reserve(Config::nBit);
    for(unsigned i = 0; i < Config::nBit; i++) {
      if(bits[i]) {
        ixv.push_back(i + 1);
        hw++;
      }
    }
    ixs = ixv;
  }

  FixedGene(size_t op_, unsigned tgt_, unsigned hw_, arma::uvec ixs_):
    op(op_), tgt(tgt_), hw(hw_), ixs(ixs_) { }

}; // class Gene

} // namespace Wrapper
