#ifndef GATE_FIXED_HPP
#define GATE_FIXED_HPP

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


template<class GateBase>
class Fixed : public GateBase {

  size_t op;
  unsigned tgt;
  Backend::Controls ixs;

  using typename GateBase::Pointer;

public:

  static Pointer getNew() {
    /* Distributions: cheap and safer in MT environment this way */
    // distribution of possible gates (except of identity)
    std::uniform_int_distribution<size_t> dOp{1, gates.size() - 1};
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    // distribution of controls
    unsigned tgt_ = dTgt(gen::rng);
    Tools::controls_distribution dCtrl{Config::nBit, Config::pControl, tgt_};
    return std::make_shared<Fixed>(
        dOp(gen::rng), tgt_, dCtrl(gen::rng));
  }

  Backend::State applyTo(const Backend::State& psi) const override {
    return psi.apply_ctrl(gates[op].op, ixs, tgt);
  }

  bool isTrivial() const override {
    return op == 0;
  }

  unsigned complexity() const override {
    return ixs.size() * ixs.size();
  }

  void invert(Pointer& self) const override {
    int dIx = gates[op].inv;
    if(dIx != 0)
      self = std::make_shared<Fixed>(op + dIx, tgt, ixs);
  }

  bool invite(Pointer& first, Pointer& second) const override {
    return first->merge(first, second, *this);
  }

  bool merge(Pointer& first, Pointer&, const Fixed& g) const override {
    // G * G = square(G) if also among our operations
    if(g.op == op && g.tgt == tgt && g.ixs == ixs && gates[op].sq != 0) {
      first = std::make_shared<Fixed>(op + gates[op].sq, tgt, ixs);
      return true;
    } else
      return false;
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

  NOINLINE Fixed(size_t op_, unsigned tgt_, std::vector<bool> ctrl):
      op(op_), tgt(tgt_), ixs(ctrl) { }

  Fixed(size_t op_, unsigned tgt_, const Backend::Controls& ixs_):
    op(op_), tgt(tgt_), ixs(ixs_) { }

}; // class Fixed

} // namespace QGA

#endif // !defined GATE_FIXED_HPP
