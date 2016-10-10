// allow only one wrapper
#ifndef QGA_WRAPPER_HPP
#define QGA_WRAPPER_HPP

#define QICLIB_DONT_USE_NLOPT
#define ARMA_DONT_USE_WRAPPER

#include "QIClib"
#include <string>
#include <cmath>
#include <sstream>

namespace Wrapper {

using Gate = arma::cx_mat22;


namespace internal {


/* Useful constants and typedefs */

using cxd = arma::cx_double;

const double pi = std::acos(-1);

const cxd i{0,1};

const double v12 = 1/std::sqrt(2);


/* Fixed gates */

Gate I {
  1, 0, 0, 1
};

Gate H {
  v12, v12, v12, -v12
};

Gate X {
  0, 1, 1, 0
};

Gate Y {
  0, -i, i, 0
};

Gate Z {
  1, 0, 0, -1
};

Gate T {
  1, 0, 0, std::exp(i*pi/4.)
};

Gate Ti {
  1, 0, 0, std::exp(-i*pi/4.)
};

Gate S {
  1, 0, 0, i
};

Gate Si {
  1, 0, 0, -i
};


/* Parametric gates */

Gate xrot(double a) {
  return { std::cos(a), i*std::sin(a), i*std::sin(a), std::cos(a) };
}

Gate yrot(double a) {
  return { std::cos(a), std::sin(a), -std::sin(a), std::cos(a) };
}

Gate zrot(double a) {
  return { std::exp(i*a), 0, 0, std::exp(-i*a) };
}

// An assymetric version of zrot
Gate phase(double a) {
  return { 1, 0, 0, std::exp(i*a) };
}


class Controls : public arma::uvec {

public:

  Controls() = default;

  Controls(const std::vector<bool>& bits) : arma::uvec(ix_vector(bits)) { }

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


} // namespace internal


class State : public arma::cx_vec {

  using Base = arma::cx_vec;
  using Controls = internal::Controls;

public:

  State(const Base& base) : Base(base) { }

  State(Base&& base) : Base(std::move(base)) { }

  // initializes in a basis state
  State(size_t index = 0) : Base(dim(), arma::fill::zeros) {
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

  static internal::cxd overlap(const State& lhs, const State& rhs) {
    return arma::cdot(lhs.rep(), rhs.rep());
  }

  template<class Gene, class... Args>
  void apply(const Gene& g, Args... args) {
    g.applyTo(*this, args...);
  }

  void apply_ctrl(const Gate& mat, const Controls& ixs, unsigned tgt) {
    *this = {qic::apply_ctrl(rep(), mat, ixs, {tgt + 1})};
  }

private:

  const Base& rep() const {
    return static_cast<const Base&>(*this);
  }

  static arma::uword dim() {
    return arma::uword(1) << Config::nBit;
  }

}; // class State

} // namespace Wrapper

#endif // !defined QGA_WRAPPER_HPP
