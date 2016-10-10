// allow only one problem
#ifndef QGA_PROBLEM_HPP
#define QGA_PROBLEM_HPP

#include "FixedGene.hpp"

namespace Wrapper {

Wrapper::State out{3};


class Gene : public QGA::GeneBase<Gene, FixedGene> {

public:

  static std::shared_ptr<Gene> getNew() {
    return Wrapper::FixedGene<Gene>::getNew();
  }

}; // class Gene


class Candidate : public QGA::CandidateBase<Candidate, Gene> {

  using Base = QGA::CandidateBase<Candidate, Gene>;

public:

  using Base::Base;

  double error() const {
    return 1 - std::abs(State::overlap(out, sim()));
  }

  std::string dump(const std::ostream& ex) const {
    std::ostringstream os{};
    os.flags(ex.flags());
    os.precision(ex.precision());
    os << sim() << '\n';
    return os.str();
  }

private:

  State sim() const {
    State psi{0};
    for(const auto& g : gt)
      psi.apply(*g);
    return psi;
  }

}; // class Candidate

} // namespace Wrapper

#endif // !defined QGA_PROBLEM_HPP
