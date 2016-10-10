// allow only one problem
#ifndef QGA_PROBLEM_HPP
#define QGA_PROBLEM_HPP

#include "XYZGene.hpp"

namespace Wrapper {

class Gene : public QGA::GeneBase<Gene, XYZGene> {

public:

  static std::shared_ptr<Gene> getNew() {
    return Wrapper::XYZGene<Gene>::getNew();
  }

}; // class Gene


class Candidate : public QGA::CandidateBase<Candidate, Gene> {

  using Base = QGA::CandidateBase<Candidate, Gene>;

public:

  using Base::Base;

  double error() const {
    if(gt.size() > 1000)
      return INFINITY;
    double error{0};
    unsigned dim = 1 << Config::nBit;
    State psi{};
    for(unsigned i = 0; i < dim; i++) {
      psi.reset(i);
      State out = State::fourier(psi);
      error += std::max(1 - std::real(State::overlap(out, sim(psi))), 0.0);
    }
    error /= dim;
    return error < 1E-8 ? 0 : error;
  }

  /*friend std::ostream& operator<< (std::ostream& os, const Candidate& c) {
    os << (Base&)c;
    double phase = 0;
    for(auto& g : c.gt)
      phase += g.phase();
    os << "φ " << phase / internal::pi << "π";
    return os;
  }*/

  std::string dump(const std::ostream& ex) const {
    std::ostringstream os{};
    os.flags(ex.flags());
    os.precision(ex.precision());
    os << '\n';
    unsigned dim = 1 << Config::nBit;
    State psi{};
    for(unsigned i = 0; i < dim; i++) {
      psi.reset(0);
      for(auto& p : sim(psi))
        os << std::abs(p)*std::sqrt(dim) << "/√" << dim << "∠"
          << std::showpos << std::arg(p)/internal::pi << "π " << std::noshowpos;
      os << '\n';
    }
    return os.str();
  }

private:

  State sim(State& psi) const {
    for(const auto& g : gt)
      psi.apply(*g);
    return psi;
  }

}; // class Candidate


} // namespace Wrapper

#endif // !defined QGA_PROBLEM_HPP
