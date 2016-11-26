#include "QGA_full.hpp"

/* Implementation of routines from QGA_bits/Gates.hpp and QGA_bits/Tools.hpp */

namespace QGA {
namespace Gates {
namespace func {

using QGA::Const::i;

/* Parametric gates */

Backend::Gate xrot(double a) {
  return {
    std::cos(a/2.0),   i*std::sin(a/2.0),
    i*std::sin(a/2.0), std::cos(a/2.0)
  };
}

Backend::Gate yrot(double a) {
  return {
    std::cos(a/2.0), -std::sin(a/2.0),
    std::sin(a/2.0), std::cos(a/2.0)
  };
}

Backend::Gate zrot(double a) {
  return {
    std::exp(i*a/2.0), 0,
    0, std::exp(-i*a/2.0)
  };
}

// Det -1 version of yrot (covers Hadamard)
// DO NOT USE in Gates::Param: does not represent a 1-parametric group!
/*Backend::Gate rrot(double a) {
  return {
    std::cos(a/2.0), std::sin(a/2.0),
    std::sin(a/2.0), -std::cos(a/2.0)
  };
}*/

// An assymetric version of zrot
Backend::Gate phase(double a) {
  return {
    1, 0,
    0, std::exp(i*a)
  };
}

} // namespace internal
} // namespace Gates


/* Convert a floating-point number to a rational approximation. This is done
 * by finding a continued fraction expression, trimming it at a random point
 * with probability proportional to the magnitude of the corresponding term,
 * and converting back. If the number is precisely rational or almost
 * rational, almost-infinite terms are capped so it can still be trimmed
 * earlier to an even shorter rational (just with a small probability). */

double rationalize(double x) {
  double a = std::abs(x);
  constexpr unsigned N = 8;
  double coeffs[N];
  unsigned t;
  for(t = 0; t < N; t++) {
    coeffs[t] = std::floor(a);
    if(coeffs[t] > 100) {
      coeffs[t++] = 100;
      break;
    }
    a = 1/(a - coeffs[t]);
  }
  std::discrete_distribution<unsigned> dStop(&coeffs[1], &coeffs[t]);
  unsigned cut = dStop(gen::rng) + 1;
  if(cut == t)
    return x;
  a = coeffs[--cut];
  while(cut > 0)
    a = coeffs[--cut] + 1/a;
  return x < 0 ? -a : a;
}

/* The same as above for angles: the variable is supposed to be 2π-periodical
 * and is replaced by a rational approximant multiple of π between -π and +π
 * (inclusive on right). */

double rationalize_angle(double a) {
  double b = a / QGA::Const::pi / 2.0 + 0.5;
  b = rationalize(b - std::floor(b));
  if(b == 0)
    b = 1;
  return (b - 0.5) * QGA::Const::pi * 2.0;
}

} // namespace QGA
