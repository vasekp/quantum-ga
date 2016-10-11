// allow only one problem
#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include "../XYZGene.hpp"

using QGA::Backend::State;


class Gene : public QGA::GeneBase<Gene, QGA::XYZGene> {

public:

  static std::shared_ptr<Gene> getNew() {
    return QGA::XYZGene<Gene>::getNew();
  }

}; // class Gene


class Candidate : public QGA::CandidateBase<Candidate, Gene> {

  using Base = QGA::CandidateBase<Candidate, Gene>;

public:

  using Base::Base;

  double error() const {
    if(gt.size() > 1000)
      return INFINITY;
    std::complex<double> avg_overlap{0};
    unsigned dim = 1 << Config::nBit;
    State psi{};
    for(unsigned i = 0; i < dim; i++) {
      psi.reset(i);
      State out = State::fourier(psi);
      avg_overlap += State::overlap(out, sim(psi));
    }
    avg_overlap /= dim;
    double error = std::max(1.0 - std::abs(avg_overlap), 0.0);
    return error < 1E-8 ? 0 : error;
  }

  std::string dump(const std::ostream& ex) const {
    std::ostringstream os{};
    os.flags(ex.flags());
    os.precision(ex.precision());
    os << '\n';
    unsigned dim = 1 << Config::nBit;
    State psi{};
    for(unsigned i = 0; i < dim; i++) {
      psi.reset(i);
      State out = sim(psi);
      for(unsigned j = 0; j < dim; j++)
        os << std::abs(out[j])*std::sqrt(dim) << "/√" << dim << "∠"
          << std::showpos << std::arg(out[j]) / QGA::Const::pi << "π "
          << std::noshowpos;
      os << '\n';
    }
    return os.str();
  }

private:

  State sim(const State& psi) const {
    State ret{psi};
    for(const auto& g : gt)
      ret = ret.apply(*g);
    return ret;
  }

}; // class Candidate

#endif // !defined PROBLEM_HPP
