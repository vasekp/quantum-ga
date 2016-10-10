// allow only one wrapper
#ifndef QGA_WRAPPER_HPP
#define QGA_WRAPPER_HPP

#define QICLIB_DONT_USE_NLOPT
#define ARMA_DONT_USE_WRAPPER

#include "qpp.h"
#include <unsupported/Eigen/FFT>
#include <string>
#include <cmath>
#include <sstream>

namespace Wrapper {

using Gate = qpp::cmat;


namespace internal {

/* Useful constants and typedefs */

using cxd = std::complex<double>;
const cxd i{0,1};
const double pi = std::acos(-1);

const Gate I = qpp::gt.Id2;
const Gate H = qpp::gt.H;
const Gate X = qpp::gt.X;
const Gate Y = qpp::gt.Y;
const Gate Z = qpp::gt.Z;
const Gate T = qpp::gt.T;
const Gate Ti = qpp::gt.T.conjugate();
const Gate S = qpp::gt.S;
const Gate Si = qpp::gt.S.conjugate();


class Controls : public std::vector<qpp::idx> {

  using Base = std::vector<qpp::idx>;

public:

  Controls() = default;

  Controls(const std::vector<bool>& bits): Base() {
    for(unsigned i = 0; i < Config::nBit; i++)
      if(bits[i])
        Base::push_back(i);
  }

  std::vector<qpp::idx> as_vector() const {
    return static_cast<const Base&>(*this);
  }

}; // class Controls


} // namespace internal


class State : public qpp::ket {

  using Base = qpp::ket;
  using Controls = internal::Controls;

public:

  State(const Base& base) : Base(base) { }

  State(Base&& base) : Base(std::move(base)) { }

  // initializes in a basis state
  State(size_t index = 0) :
    Base(qpp::mket(qpp::n2multiidx(index, dims()))) { }

  // resets in a basis state
  void reset(size_t index) {
    *this = qpp::mket(qpp::n2multiidx(index, dims()));
  }

  static State fourier(const State& in) {
    Eigen::VectorXcd ret{in.size()};
    Eigen::FFT<double> fft;
    ret.col(0) = fft.fwd(in.col(0));
    return {ret};
  }

  static internal::cxd overlap(const State& lhs, const State& rhs) {
    return rhs.rep().dot(lhs.rep());
  }

  template<class Gene, class... Args>
  void apply(const Gene& g, Args... args) {
    g.applyTo(*this, args...);
  }

  void apply_ctrl(const Gate& mat, const Controls& ixs, unsigned tgt) {
    *this = {qpp::applyCTRL(rep(), mat, ixs, {tgt})};
  }

  friend std::ostream& operator<< (std::ostream& os, const State& state) {
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols,
        " ", " "); // row, col separators
    return os << state.format(fmt);
  }

private:

  const Base& rep() const {
    return static_cast<const Base&>(*this);
  }

  std::vector<qpp::idx> dims() {
    return std::vector<qpp::idx>(Config::nBit, 2);
  }

}; // class State


} // namespace Wrapper

#endif // !defined QGA_WRAPPER_HPP
