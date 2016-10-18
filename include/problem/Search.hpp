// allow only one problem
#ifndef QGA_PROBLEM_HPP
#define QGA_PROBLEM_HPP

namespace {

using QGA::Backend::State;

/* An extension of QGA::GateBase allowing us to count oracle calls and pass
 * an additional parameter to applyTo(). */

template<class GateBase, class... Gates>
class NewBase : public QGA::GateBase<GateBase, Gates...> {

  using QGA::GateBase<GateBase, Gates...>::applyTo;

public:

  virtual State applyTo(const State& psi, unsigned) const {
    return this->applyTo(psi);
  }

};


/* The oracle gate. */

struct Oracle {

template<class GateBase>
class Inner : public GateBase {

  using typename GateBase::Pointer;
  using typename GateBase::Counter;

  bool odd;  // parity of the power

public:

  static Pointer getNew() {
    return std::make_shared<Inner>();
  }

  State applyTo(const State&) const override {
    throw std::logic_error("Oracle::applyTo called without mark");
  }

  State applyTo(const State& psi, unsigned mark) const override {
    State ret{psi};
    if(odd)
      ret[mark] = -ret[mark];
    return ret;
  }

  bool isTrivial() const override {
    // oracle^(2k) = oracle^0 = identity
    return !odd;
  }

  unsigned complexity() const override {
    return 1;
  }

  void hit(Counter& c) const {
    c.hit(this);
  }

  Pointer invite(const Pointer& first) const override {
    return first->merge(*this);
  }

  Pointer merge(const Inner& g) const override {
    // oracle * oracle = oracle^2 → true ^ true = false
    return std::make_shared<Inner>(odd ^ g.odd);
  }

  std::ostream& write(std::ostream& os) const override {
    return os << (odd ? "Oracle" : "[Id]");
  }

  static Pointer read(const std::string& s) {
    std::regex re{"\\[Id\\]|(Oracle)"};
    std::smatch m{};
    if(!std::regex_match(s, m, re))
      return {};
    return std::make_shared<Inner>(m[1].matched);
  }

  Inner(bool odd_ = true): odd(odd_) { }

}; // class Oracle::Inner

template<class GateBase>
using Template = Inner<GateBase>;

}; // struct Oracle


/* Our Gene type will randomly choose between uncontrolled X/Y/Z, controlled
 * Phase, and Oracle and will support applyTo(const State&, unsigned) */

using Gene = QGA::CustomGene<NewBase,
        QGA::Gates::X<>,
        QGA::Gates::CPhase,
        Oracle>;


class Candidate : public QGA::CandidateBase<Candidate, Gene> {

  using Base = QGA::CandidateBase<Candidate, Gene>;

public:

  using Base::Base;

  double error() const {
    if(gt.size() > 1000)
      return INFINITY;
    double error{0};
    unsigned dim = 1 << Config::nBit;
    State psi{0};
    for(unsigned mark = 0; mark < dim; mark++) {
      State out{mark};
      error += std::max(1 -
          std::pow(std::abs(State::overlap(out, sim(psi, mark))), 2), 0.0);
    }
    return error / dim;
  }

  std::string dump(const std::ostream& ex) const {
    std::ostringstream os{};
    os.flags(ex.flags());
    os.precision(ex.precision());
    os << '\n';
    unsigned dim = 1 << Config::nBit;
    State psi{0};
    for(unsigned mark = 0; mark < dim; mark++) {
      os << mark << ": ";
      os << sim(psi, mark);
    }
    return os.str();
  }

private:

  State sim(const State& psi, unsigned mark) const {
    State ret{psi};
    for(const auto& g : gt)
      ret = ret.apply(g, mark);
    return ret;
  }

}; // class Candidate

} // anonymous namespace

#endif // !defined QGA_PROBLEM_HPP
