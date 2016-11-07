// allow only one problem
#ifndef QGA_PROBLEM_HPP
#define QGA_PROBLEM_HPP

namespace {

using QGA::Backend::State;


/* The Context to be used in Gene below, holding a mark for the Oracle */
struct Context {
  unsigned mark;
};


/* The oracle gate. */

struct Oracle {

template<class GateBase>
class OracleTemp : public GateBase {

  using typename GateBase::Pointer;
  using typename GateBase::Context;

  bool odd;  // parity of the power

public:

  OracleTemp(bool odd_ = true): odd(odd_) { }

  State applyTo(const State& psi, const Context* pMark) const override {
    unsigned mark = pMark->mark;
    State ret{psi};
    if(odd)
      ret[mark] = -ret[mark];
    return ret;
  }

  bool isTrivial() const override {
    // oracle^(2k) = oracle^0 = identity
    return !odd;
  }

  void hit(typename GateBase::Counter& c) const {
    c.hit(this);
  }

  Pointer invite(const Pointer& first) const override {
    return first->merge(*this);
  }

  Pointer merge(const OracleTemp& g) const override {
    // oracle * oracle = oracle^2 → true ^ true = false
    return std::make_shared<OracleTemp>(odd ^ g.odd);
  }

  std::ostream& write(std::ostream& os) const override {
    return os << (odd ? "Oracle" : "[Id]");
  }

  void printOn(QGA::internal::CircuitPrinter& p) const override {
    p.addBarrierGate("U_f");
  }

  static Pointer read(const std::string& s) {
    std::regex re{"\\[Id\\]|(Oracle)"};
    std::smatch m{};
    if(!std::regex_match(s, m, re))
      return {};
    return std::make_shared<OracleTemp>(m[1].matched);
  }

}; // class Oracle::OracleTemp<GateBase>

template<class GateBase>
using Template = OracleTemp<GateBase>;

}; // struct Oracle


using Gene = typename QGA::Gene<
                Oracle,
                QGA::Gates::X,
                QGA::Gates::CPhase
              >::WithContext<Context>;


class Candidate : public QGA::CandidateBase<Candidate, Gene, double, double> {

  using Base = QGA::CandidateBase<Candidate, Gene, double, double>;

public:

  using Base::Base;

  Base::FitnessMain fitness_main() const {
    if(genotype().size() > 1000)
      return {INFINITY, INFINITY};
    double errTotal = 0, errMax = 0;
    unsigned dim = 1 << Config::nBit;
    State psi{0};
    for(unsigned mark = 0; mark < dim; mark++) {
      State out{mark};
      double error = std::max(1 -
          std::pow(std::abs(State::overlap(out, sim(psi, mark))), 2), 0.0);
      errTotal += error;
      if(error > errMax)
        errMax = error;
    }
    return {
      this->trimError(errTotal / dim),
      this->trimError(errMax)
    };
  }

  std::ostream& print_full(std::ostream& os) const {
    unsigned dim = 1 << Config::nBit;
    State psi{0};
    os << '\n';
    for(unsigned mark = 0; mark < dim; mark++) {
      os << mark << ": ";
      os << sim(psi, mark);
    }
    return os;
  }

private:

  State sim(const State& psi, unsigned mark) const {
    State ret{psi};
    Context c{mark};
    for(const auto& g : genotype())
      ret = ret.apply(g, &c);
    return ret;
  }

}; // class Candidate

} // anonymous namespace

#endif // !defined QGA_PROBLEM_HPP