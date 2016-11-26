/* Common definitions */
namespace QGA {
namespace Gates {
namespace func {

using QGA::Const::i;

/* Parametric gates */

Backend::Gate xrot(double);

Backend::Gate yrot(double);

Backend::Gate zrot(double);

// An assymetric version of zrot
Backend::Gate phase(double);

} // namespace internal
} // namespace Gates
} // namespace QGA

#include "gates/Fixed.hpp"
#include "gates/XYZ.hpp"
#include "gates/CPhase.hpp"
#include "gates/CNOT.hpp"
#include "gates/SU2.hpp"
#include "gates/SWAP.hpp"
