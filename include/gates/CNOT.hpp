namespace QGA {

namespace Gates {

using Tools::Controls;

namespace internal {

template<class GateBase, Controls cc>
class CNOT : public GateBase {

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
    return std::make_shared<CNOT>(tgt, dCtrl(gen::rng));
  }

  Backend::State applyTo(const Backend::State& psi) const override {
    return odd ? psi.apply_ctrl(Backend::X, ixs, tgt) : psi;
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

  bool merge(Pointer& first, Pointer&, const CNOT& g) const override {
    if(g.tgt == tgt && g.ixs == ixs) {
      first = std::make_shared<CNOT>(tgt, ixs, odd ^ g.odd);
      return true;
    } else
      return false;
  }

  std::ostream& write(std::ostream& os) const override {
    if(!odd)
      return os << "[Id]";
    os << "NOT" << tgt + 1;
    if(ixs.size()) {
      os << '[';
      for(auto ctrl : ixs.as_vector())
        os << ctrl + 1;
      os << ']';
    }
    return os;
  }

  CNOT(unsigned tgt_, std::vector<bool> ctrl):
      tgt(tgt_), ixs(ctrl), odd(true) { }

  CNOT(unsigned tgt_, const Backend::Controls& ixs_, bool odd_ = true):
      tgt(tgt_), ixs(ixs_), odd(odd_) { }

}; // class CNOT

} // namespace internal


template<class GateBase>
using CNOT = internal::CNOT<GateBase, Controls::ONE>;

template<class GateBase>
using CkNOT = internal::CNOT<GateBase, Controls::LEAST1>;

} // namespace Gates

} // namespace QGA
