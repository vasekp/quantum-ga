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


/* An enum for possible settings for controls_distribution below */

enum class Controls {
  NONE,
  ONE,
  LEAST1,
  ANY
};

/* A distribution generating bit strings of length nBit where the probability
 * of 1 in each position is given by pTrue. The bit at position nSkip is left
 * off. */

template<Controls cc>
class controls_distribution {

  const unsigned nBit;
  const double pTrue;
  const unsigned iSkip;

public:

  controls_distribution(unsigned nBit_, unsigned iSkip_, double pTrue_):
    nBit(nBit_), pTrue(pTrue_), iSkip(iSkip_) {
#ifdef DEBUG
      if((cc == Controls::ONE || cc == Controls::LEAST1) && nBit <= 1)
        throw std::logic_error("nBit < 2 in controls_distribution<≥1>!");
#endif
    }

  template<class URNG>
  std::vector<bool> operator() (URNG& rng) {
    std::vector<bool> bits(nBit, false);

    if(cc == Controls::ANY || cc == Controls::LEAST1) {
      std::bernoulli_distribution dist(pTrue);
      for(unsigned i = 0; i < nBit; i++) {
        if(i == iSkip)
          continue;
        bits[i] = dist(rng);
      }
    }

    if(cc == Controls::ONE || cc == Controls::LEAST1) {
      std::uniform_int_distribution<> dist(0, nBit - 2);
      unsigned res = dist(rng);
      bits[res + (res >= iSkip)] = true;
    }

    return bits;
  }

}; // class controls_distribution<Controls>

} // namespace Tools

} // namespace QGA
