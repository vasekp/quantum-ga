// allow only one backend
#ifndef QGA_BACKEND_HPP
#define QGA_BACKEND_HPP

#define QICLIB_DONT_USE_NLOPT
#define ARMA_DONT_USE_WRAPPER

#include "QIClib"

namespace QGA {

namespace Backend {

using Gate = arma::cx_mat22;

using QGA::Const::i;
using QGA::Const::pi;
using QGA::Const::v12;

/* Fixed gates */

const Gate I { 1, 0, 0, 1 };
const Gate H { v12, v12, v12, -v12 };
const Gate X { 0, 1, 1, 0 };
const Gate Y { 0, -i, i, 0 };
const Gate Z { 1, 0, 0, -1 };
const Gate T { 1, 0, 0, std::exp(i*pi/4.) };
const Gate Ti { 1, 0, 0, std::exp(-i*pi/4.) };
const Gate S { 1, 0, 0, i };
const Gate Si { 1, 0, 0, -i };

/* Parametric gates */

Gate xrot(double a) {
  return {
    std::cos(a/2.0),   i*std::sin(a/2.0),
    i*std::sin(a/2.0), std::cos(a/2.0)
  };
}

Gate yrot(double a) {
  return {
    std::cos(a/2.0), -std::sin(a/2.0),
    std::sin(a/2.0), std::cos(a/2.0)
  };
}

Gate zrot(double a) {
  return {
    std::exp(i*a/2.0), 0,
    0, std::exp(-i*a/2.0)
  };
}

// An assymetric version of zrot
Gate phase(double a) {
  return {
    1, 0,
    0, std::exp(i*a)
  };
}


class Controls : public arma::uvec {

public:

  Controls() = default;

  Controls(const std::vector<bool>& bits): arma::uvec(ix_vector(bits)) { }

  friend bool operator== (const Controls& lhs, const Controls& rhs) {
    return lhs.size() == rhs.size() && arma::all(lhs.rep() == rhs.rep());
  }

  std::vector<unsigned> as_vector() const {
    std::vector<unsigned> ret{};
    ret.reserve(size());
    for(auto ix : *this)
      ret.push_back(ix - 1);
    return ret;
  }

private:

  const arma::uvec& rep() const {
    return static_cast<const arma::uvec&>(*this);
  }

  static std::vector<arma::uword> ix_vector(const std::vector<bool>& bits) {
    std::vector<arma::uword> ret{};
    for(unsigned i = 0; i < Config::nBit; i++)
      if(bits[i])
        ret.push_back(i + 1);
    return ret;
  }

}; // class Controls


class State : public arma::cx_vec {

  using Base = arma::cx_vec;

public:

  State(const Base& base): Base(base) { }

  State(Base&& base): Base(std::move(base)) { }

  // initializes in a basis state
  State(size_t index = 0): Base(dim(), arma::fill::zeros) {
    this->operator[](index) = 1;
  }

  // resets in a basis state
  void reset(size_t index) {
    this->fill(0);
    this->operator[](index) = 1;
  }

  static State fourier(const State& in) {
    return {arma::fft(in.rep()) / sqrt(dim())};
  }

  static std::complex<double> overlap(const State& lhs, const State& rhs) {
    return arma::cdot(lhs.rep(), rhs.rep());
  }

  template<class Gene, class... Args>
  State apply(const Gene& g, Args... args) {
    return g->applyTo(*this, args...);
  }

  State apply_ctrl(const Gate& mat, const Controls& ixs, unsigned tgt) const {
    return {qic::apply_ctrl(rep(), mat, ixs, {tgt + 1})};
  }

  friend std::ostream& operator<< (std::ostream& os, const State& state) {
    state.st().raw_print(os);
    return os;
  }

private:

  const Base& rep() const {
    return static_cast<const Base&>(*this);
  }

  static arma::uword dim() {
    return arma::uword(1) << Config::nBit;
  }

}; // class State

} // namespace Backend

} // namespace QGA

#endif // !defined QGA_BACKEND_HPP
