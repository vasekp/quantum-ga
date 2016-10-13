#ifndef GATE_CNOT_HPP
#define GATE_CNOT_HPP

namespace QGA {

using Tools::Controls;

template<class GateBase, Controls cc>
class CNot : public GateBase {

  unsigned tgt;
  Backend::Controls ixs;
  bool odd;  // parity of the power

  using typename GateBase::Pointer;

public:

  static Pointer getNew() {
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    // distribution of controls
    unsigned tgt = dTgt(gen::rng);
    Tools::controls_distribution<cc> dCtrl{Config::nBit, tgt, Config::pControl};
    return std::make_shared<CNot>(tgt, dCtrl(gen::rng));
  }

  Backend::State applyTo(const Backend::State& psi) const override {
    return psi.apply_ctrl(Backend::X, ixs, tgt);
  }

  bool isTrivial() const override {
    // CNOT^(2k) = CNOT^0 = identity
    return !odd;
  }

  unsigned complexity() const override {
    return ixs.size() * ixs.size();
  }

  bool invite(Pointer& first, Pointer& second) const override {
    return first->merge(first, second, *this);
  }

  bool merge(Pointer& first, Pointer&, const CNot& g) const override {
    if(g.tgt == tgt && g.ixs == ixs) {
      first = std::make_shared<CNot>(tgt, ixs, odd ^ g.odd);
      return true;
    } else
      return false;
  }

  std::ostream& write(std::ostream& os) const override {
    os << "NOT" << tgt + 1;
    if(ixs.size()) {
      os << '[';
      for(auto ctrl : ixs.as_vector())
        os << ctrl + 1;
      os << ']';
    }
    return os;
  }

  CNot(unsigned tgt_, std::vector<bool> ctrl):
      tgt(tgt_), ixs(ctrl), odd(true) { }

  CNot(unsigned tgt_, const Backend::Controls& ixs_, bool odd_ = true):
      tgt(tgt_), ixs(ixs_), odd(odd_) { }

}; // class CNot

template<class GateBase>
using CNOT = CNot<GateBase, Controls::ONE>;

template<class GateBase>
using CnNOT = CNot<GateBase, Controls::ANY>;

} // namespace QGA

#endif // !defined GATE_CNOT_HPP
