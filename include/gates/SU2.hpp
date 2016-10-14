#ifndef GATE_SU2_HPP
#define GATE_SU2_HPP

namespace QGA {

namespace Gates {

using Tools::Controls;

namespace internal {

template<class GateBase, Controls cc>
class SU2 : public GateBase {

  unsigned tgt;
  double angle1;
  double angle2;
  double angle3;
  Backend::Controls ixs;
  Backend::Gate mat;

  using typename GateBase::Pointer;

public:

  static Pointer getNew() {
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    // distribution of controls
    unsigned tgt = dTgt(gen::rng);
    Tools::controls_distribution<cc> dCtrl{Config::nBit, tgt, Config::pControl};
    // distribution of angle
    std::uniform_real_distribution<> dAng{-0.5*Const::pi, 0.5*Const::pi};
    return std::make_shared<SU2>(tgt,
        dAng(gen::rng), dAng(gen::rng), dAng(gen::rng), dCtrl(gen::rng));
  }

  Backend::State applyTo(const Backend::State& psi) const override {
    return psi.apply_ctrl(mat, ixs, tgt);
  }

  bool isTrivial() const override {
    return angle2 == 0 && angle1 + angle3 == 0;
  }

  unsigned complexity() const override {
    return ixs.size() * ixs.size();
  }

  void invert(Pointer& self) const override {
    self = std::make_shared<SU2>(tgt, -angle3, -angle2, -angle1, ixs);
  }

  void mutate(Pointer& self) const override {
    std::normal_distribution<> dAng{0.0, 0.1};
    self = std::make_shared<SU2>(tgt,
        angle1 + dAng(gen::rng),
        angle2 + dAng(gen::rng),
        angle3 + dAng(gen::rng),
        ixs);
  }

  void simplify(Pointer& self) const override {
    self = std::make_shared<SU2>(tgt,
        Tools::rationalize_angle(angle1),
        Tools::rationalize_angle(angle2),
        Tools::rationalize_angle(angle3),
        ixs);
  }

  bool invite(Pointer& first, Pointer& second) const override {
    return first->merge(first, second, *this);
  }

  bool merge(Pointer& first, Pointer&, const SU2& g) const override {
    if(g.tgt == tgt && g.ixs == ixs) {
      first = std::make_shared<SU2>(tgt, ixs, static_cast<Backend::Gate>(g.mat * mat));
      return true;
    } else
      return false;
  }

  std::ostream& write(std::ostream& os) const override {
    os << "U" << tgt + 1;
    if(ixs.size()) {
      os << '[';
      for(auto ctrl : ixs.as_vector())
        os << ctrl + 1;
      os << ']';
    }
    os << '('
       << angle1 / Const::pi << "π,"
       << angle2 / Const::pi << "π,"
       << angle3 / Const::pi << "π)";
    return os;
  }

  SU2(unsigned tgt_, double angle1_, double angle2_, double angle3_,
      std::vector<bool> ctrl):
    tgt(tgt_), angle1(angle1_), angle2(angle2_), angle3(angle3_), ixs(ctrl),
    mat(Backend::zrot(angle3) * Backend::yrot(angle2) * Backend::zrot(angle1))
  { }

  SU2(unsigned tgt_, double angle1_, double angle2_, double angle3_,
      const Backend::Controls& ixs_):
    tgt(tgt_), angle1(angle1_), angle2(angle2_), angle3(angle3_), ixs(ixs_),
    mat(Backend::zrot(angle3) * Backend::yrot(angle2) * Backend::zrot(angle1))
  { }

  SU2(unsigned tgt_, const Backend::Controls& ixs_, Backend::Gate&& mat_):
    tgt(tgt_), angle1(), angle2(), angle3(), ixs(ixs_), mat(mat_)
  {
    angle2 = std::atan2(std::abs(mat(1, 0)), std::abs(mat(0, 0)));
    double sum = std::arg(mat(0, 0)),
           diff = std::arg(mat(1, 0));
    angle1 = (sum + diff) / 2.0;
    angle3 = (sum - diff) / 2.0;
  }

}; // class SU2

} // namespace internal


template<class GateBase>
using SU2 = internal::SU2<GateBase, Controls::NONE>;

} // namespace Gates

} // namespace QGA

#endif // !defined GATE_SU2_HPP
