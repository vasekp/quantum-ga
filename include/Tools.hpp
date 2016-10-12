namespace QGA {

namespace Tools {

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
 * and is replaced by a rational approximant multiple of π. */

double rationalize_angle(double a) {
  return rationalize(std::fmod(a / QGA::Const::pi, 2.0)) * QGA::Const::pi;
}


/* A distribution generating bit strings of length nBit where the probability
 * of 1 in each position is given by pTrue. The bit at position nSkip is left
 * off. The provided uniform RNG is invoked only the minimum necessary number
 * of times. */
class controls_distribution {

  const unsigned nBit;
  const double pTrue;
  const unsigned iSkip;

  // entropy decrease per true result
  const double dSTrue;
  // entropy decrease per false result
  const double dSFalse;
  // entropy minimum left for reliable generation
  const double Smin;

public:

  controls_distribution(unsigned nBit_, double pTrue_,
      unsigned iSkip_ = (unsigned)(~0)):
    nBit(nBit_), pTrue(pTrue_), iSkip(iSkip_),
    dSTrue(-std::log(pTrue) / std::log(2.0)),
    dSFalse(-std::log(1 - pTrue) / std::log(2.0)),
    Smin(std::max(dSTrue, dSFalse)) { }

  template<class URNG>
  std::vector<bool> operator() (URNG& rng) {
    std::vector<bool> bits(nBit, false);
    std::uniform_real_distribution<> dist;
    double c = dist(rng);
    // available entropy per call to dist
    constexpr double S0 =
      std::numeric_limits<typename URNG::result_type>::digits;
    // current entropy
    double S = S0;

    for(unsigned i = 0; i < nBit; i++) {
      if(i == iSkip)
        continue;
      if(c < pTrue) {
        bits[i] = true;
        c /= pTrue;
        S -= dSTrue;
      } else {
        c = (c - pTrue)/(1 - pTrue);
        S -= dSFalse;
      }
      if(S < Smin) {
        c = dist(rng);
        S = S0; // NB can still be < Smin but there's not much more we can do
      }
    }
    return bits;
  }

}; // class controls_distribution

} // namespace Tools

} // namespace QGA
