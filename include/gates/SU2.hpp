namespace QGA {

namespace Gates {

namespace internal {

template<Controls cc>
struct SU2 {

template<class GateBase>
class SU2Temp : public GateBase {

  unsigned tgt;
  double angle1;
  double angle2;
  double angle3;
  Backend::Controls ixs;
  Backend::Gate mat;

  using typename GateBase::Pointer;
  using Ctx = typename GateBase::Context;

public:

  // construct a random gate
  SU2Temp():
    tgt(std::uniform_int_distribution<unsigned>{0, Config::nBit - 1}
        (gen::rng)),
    angle1(angle_distribution<>{}(gen::rng)),
    angle2(angle_distribution<>{}(gen::rng)),
    angle3(angle_distribution<>{}(gen::rng)),
    ixs(controls_distribution<cc>{Config::nBit, tgt, Config::pControl}
        (gen::rng)),
    mat(func::zrot(angle3) * func::yrot(angle2) * func::zrot(angle1))
  { }

  // construct using parameters
  SU2Temp(unsigned tgt_, double angle1_, double angle2_, double angle3_,
      const Backend::Controls& ixs_):
    tgt(tgt_), angle1(angle1_), angle2(angle2_), angle3(angle3_), ixs(ixs_),
    mat(func::zrot(angle3) * func::yrot(angle2) * func::zrot(angle1))
  { }

  // construct from a product matrix
  SU2Temp(unsigned tgt_, const Backend::Controls& ixs_, Backend::Gate&& mat_):
    tgt(tgt_), angle1(), angle2(), angle3(), ixs(ixs_), mat(mat_)
  {
    angle2 = std::atan2(std::abs(mat(1, 0)), std::abs(mat(0, 0)));
    double sum = std::arg(mat(0, 0)),
           diff = std::arg(mat(1, 0));
    angle1 = (sum + diff) / 2.0;
    angle3 = (sum - diff) / 2.0;
  }

  Backend::State applyTo(const Backend::State& psi, const Ctx*) const override {
    return psi.apply_ctrl(mat, ixs, tgt);
  }

  bool isTrivial() const override {
    return angle2 == 0 && angle1 + angle3 == 0;
  }

  unsigned controls() const override {
    return ixs.size();
  }

  Pointer invert(const Pointer&) const override {
    return std::make_shared<SU2Temp>(tgt, -angle3, -angle2, -angle1, ixs);
  }

  Pointer mutate(const Pointer&) const override {
    angle_distribution<true> dAng{};
    return std::make_shared<SU2Temp>(tgt,
        angle1 + dAng(gen::rng),
        angle2 + dAng(gen::rng),
        angle3 + dAng(gen::rng),
        ixs);
  }

  Pointer simplify(const Pointer&) const override {
    return std::make_shared<SU2Temp>(tgt,
        rationalize_angle(angle1),
        rationalize_angle(angle2),
        rationalize_angle(angle3),
        ixs);
  }

  void hit(typename GateBase::Counter& c) const {
    c.hit(this);
  }

  Pointer invite(const Pointer& first) const override {
    return first->merge(*this);
  }

  Pointer merge(const SU2Temp& g) const override {
    if(g.tgt == tgt && g.ixs == ixs)
      return std::make_shared<SU2Temp>(tgt, ixs,
          static_cast<Backend::Gate>(g.mat * mat));
    else
      return {};
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

  static Pointer read(const std::string& s) {
    std::string reS{};
    std::regex re{"U(\\d)(\\[(\\d+)\\])?"
      "\\((-?[0-9.]+)(?:π)?,(-?[0-9.]+)(?:π)?,(-?[0-9.]+)(?:π)?\\)"};
    std::smatch m{};
    if(!std::regex_match(s, m, re))
      return {};
    unsigned tgt = m[1].str()[0] - '1';
    if(tgt >= Config::nBit)
      return {};
    std::vector<bool> ctrl(Config::nBit, false);
    if(m[2].matched)
      for(const char& c : m[3].str()) {
        size_t pos = c - '1';
        if(pos >= 0 && pos < Config::nBit && pos != tgt)
          ctrl[pos] = true;
      }
    double angle1 = std::stod(m[4].str()) * Const::pi,
           angle2 = std::stod(m[5].str()) * Const::pi,
           angle3 = std::stod(m[6].str()) * Const::pi;
    return std::make_shared<SU2Temp>(tgt,
        angle1, angle2, angle3,
        Backend::Controls{ctrl});
  }

}; // class SU2<Controls>::SU2Temp<GateBase>

template<class GateBase>
using Template = SU2Temp<GateBase>;

template<Controls cc_>
using WithControls = SU2<cc_>;

}; // struct SU2<Controls>

} // namespace internal

using SU2 = internal::SU2<Controls::NONE>;

} // namespace Gates

} // namespace QGA
