// allow only one problem
#ifndef QGA_PROBLEM_HPP
#define QGA_PROBLEM_HPP

#include "../gates/XPhase.hpp"

using QGA::Backend::State;

/* An extension of QGA::GateBase allowing us to count oracle calls and pass
 * an additional parameter to applyTo(). */

template<class GateBase, template<class> class... Gates>
class NewBase : public QGA::GateBase<GateBase, Gates...> {

  using QGA::GateBase<GateBase, Gates...>::applyTo;

public:

  virtual unsigned calls() const {
    return 0;
  }

  virtual State applyTo(const State& psi, unsigned) const {
    return this->applyTo(psi);
  }

};


/* The oracle gene template. */

template<class GateBase>
class Oracle : public GateBase {

  using typename GateBase::Pointer;
  bool odd;  // parity of the power

public:

  static Pointer getNew() {
    return std::make_shared<Oracle>();
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

  unsigned calls() const override {
    return 1;
  }

  bool invite(Pointer& first, Pointer& second) const override {
    return first->merge(first, second, *this);
  }

  bool merge(Pointer& first, Pointer&, const Oracle& g) const override {
    // oracle * oracle = oracle^2 → true ^ true = false
    first = std::make_shared<Oracle>(odd ^ g.odd);
    return true;
  }

  std::ostream& write(std::ostream& os) const override {
    return os << (odd ? "Oracle" : "[Id]");
  }

  Oracle(bool odd_ = true): odd(odd_) { }

}; // class Oracle<GateBase>


/* Our Gene type will randomly choose between XPhase and Oracle and will support
 * Gene::calls(). */

using Gene = QGA::CustomGene<NewBase, QGA::X, QGA::CPhase, Oracle>;


struct Fitness {

  double error;
  size_t length;
  size_t ocalls;

  friend std::ostream& operator<< (std::ostream& os, const Fitness& f) {
    return os << '{'
       << f.error << ','
       << f.length << ','
       << f.ocalls << '}';
  }

  friend NOINLINE bool operator< (const Fitness& a, const Fitness& b) {
    return a.error < b.error || (a.error == b.error && a.ocalls < b.ocalls);
  }

  friend NOINLINE bool operator<< (const Fitness& a, const Fitness& b) {
    return a.error <= b.error
        && a.length <= b.length
        && a.ocalls <= b.ocalls
        && !(a == b);
  }

  friend bool operator== (const Fitness& a, const Fitness& b) {
    return a.error == b.error
        && a.length == b.length
        && a.ocalls == b.ocalls;
  }

}; // struct Fitness


class Candidate : public QGA::CandidateBase<Candidate, Gene> {

  using Base = QGA::CandidateBase<Candidate, Gene>;

public:

  using Base::Base;

  NOINLINE Fitness fitness() const {
    size_t ocalls = 0;
    for(const auto& g : gt)
      ocalls += g->calls();
    QGA::counter.hit();
    return {trimError(error()), gt.size(), ocalls};
  }

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

#endif // !defined QGA_PROBLEM_HPP
