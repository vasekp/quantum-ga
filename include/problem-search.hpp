// allow only one problem
#ifndef QGA_PROBLEM_HPP
#define QGA_PROBLEM_HPP

#include "XYZGene.hpp"

namespace Wrapper {


template<class GeneBase>
class Oracle : public GeneBase {

  using SP = std::shared_ptr<GeneBase>;

public:

  Oracle() { }

  static SP getNew() {
    return std::make_shared<Oracle>();
  }

  void applyTo(State&) const override {
    throw std::logic_error("Oracle::applyTo called without mark");
  }

  void applyTo(State& psi, unsigned mark) const override {
    psi[mark] = -psi[mark];
  }

  unsigned complexity() const override {
    return 1;
  }

  unsigned calls() const override {
    return 1;
  }

  SP invite(const SP& g) const override {
    return g->visit(g, *this);
  }

  std::ostream& write(std::ostream& os) const override {
    return os << "Oracle";
  }

}; // class Oracle<GeneBase>


class Gene : public QGA::GeneBase<Gene, Oracle, XYZGene> {

public:

  using QGA::GeneBase<Gene, Oracle, XYZGene>::applyTo;

  virtual unsigned calls() const {
    return 0;
  }

  virtual void applyTo(State& psi, unsigned) const {
    return applyTo(psi);
  }

  static std::shared_ptr<Gene> getNew() {
    std::bernoulli_distribution dOracle{0.1};
    if(dOracle(gen::rng))
      return Oracle<Gene>::getNew();
    else
      return XYZGene<Gene>::getNew();
  }

}; // class GeneBase<Derived...>


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
    return {error(), gt.size(), ocalls};
  }

  double error() const {
    if(gt.size() > 1000)
      return INFINITY;
    double error{0};
    unsigned dim = 1 << Config::nBit;
    State psi{}, out{};
    for(unsigned mark = 0; mark < dim; mark++) {
      psi.reset(0);
      out.reset(mark);
      error += std::max(1 -
          std::pow(std::abs(State::overlap(out, sim(psi, mark))), 2), 0.0);
    }
    error /= dim;
    return (unsigned long)(error * (1UL<<24)) / (double)(1UL<<24);
  }

  std::string dump(const std::ostream& ex) const {
    std::ostringstream os{};
    os.flags(ex.flags());
    os.precision(ex.precision());
    os << '\n';
    unsigned dim = 1 << Config::nBit;
    State psi{};
    for(unsigned mark = 0; mark < dim; mark++) {
      os << mark << ": ";
      psi.reset(0);
      os << sim(psi, mark) << '\n';
    }
    return os.str();
  }

private:

  State sim(State& psi, unsigned mark) const {
    for(const auto& g : gt)
      psi.apply(*g, mark);
    return psi;
  }

}; // class Candidate


} // namespace Wrapper

#endif // !defined QGA_PROBLEM_HPP
