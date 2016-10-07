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

using State = arma::cx_vec;
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

} // namespace internal

} // namespace Wrapper

#endif // !defined QGA_WRAPPER_HPP
