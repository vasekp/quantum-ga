/* Common definitions */
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
Backend::Gate rrot(double a) {
  return {
    std::cos(a/2.0), std::sin(a/2.0),
    std::sin(a/2.0), -std::cos(a/2.0)
  };
}

// An assymetric version of zrot
Backend::Gate phase(double a) {
  return {
    1, 0,
    0, std::exp(i*a)
  };
}

} // namespace internal
} // namespace Gates
} // namespace QGA

#include "gates/Fixed.hpp"
#include "gates/XYZ.hpp"
#include "gates/CPhase.hpp"
#include "gates/CNOT.hpp"
#include "gates/SU2.hpp"
#include "gates/SWAP.hpp"
