namespace QGA {

/* Forward declarations */
namespace Backend {
  class State;
}

template<class, template<class> class, template<class> class...>
class Visitors;
/* End forward declarations */


/* The base class for all genes. Defines methods derived classes have to
 * implement, and provides default (no-op) definition for some of them. */
template<class Gene,
  template<class> class... Derived>
class GeneBase : public Visitors<Gene, Derived...> {

  using SP = std::shared_ptr<Gene>;

public:

  // apply this gene to a state vector
  virtual Backend::State applyTo(const Backend::State&) const = 0;

  // return an arbitrary notion of complexity of this operation (accumulative)
  virtual unsigned complexity() const = 0;

  /* In the following functions, self is a std::shared_ptr (SP) to this (so it
   * holds that self.get() == this). This means the same information is passed
   * to the function twice, but the SP could not be retrieved from this alone.
   * If the function did not modify the state it is required to return self,
   * otherwise a shared pointer to a new instance. */

  virtual SP invert(const SP& self) {
    return self;
  }

  virtual SP mutate(const SP& self) {
    return self;
  }

  virtual SP simplify(const SP& self) {
    return self;
  }

  /* Merge needs to be implemented using double dispatch, because only
   * genes of the same class can typically be merged (albeit this is not a
   * restriction) but all we know at compile time is a pointer to the base
   * class.
   *
   * See: http://www.oodesign.com/visitor-pattern.html */

  SP merge(const SP& self, const SP& other) {
    return other.get()->invite(self);
  }

  friend std::ostream& operator<< (std::ostream& os, const GeneBase& g) {
    return g.write(os);
  }

protected:

  virtual SP invite(const SP&) const = 0;

  virtual std::ostream& write(std::ostream&) const = 0;

  /* Convert a floating-point number to a rational approximation. This is done
   * by finding a continued fraction expression, trimming it at a random point
   * with probability proportional to the magnitude of the corresponding term,
   * and converting back. If the number is precisely rational or almost
   * rational, almost-infinite terms are capped so it can still be trimmed
   * earlier to an even shorter rational (just with a small probability). */
  static double rationalize(double x) {
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

}; // virtual class GeneBase<Gene, Derived...>


/* The purpose of this class is to inject a virtual method for calling a
 * particular gene class. Note that for the visitor design pattern, we need
 * one method like this for each possible visitee.
 *
 * This is for the purposes of merge(). The default implementation returns
 * self, meaning the two genes could not be merged. Subclasses can override
 * SP visit(const SP&, const X<Gene>&) for some X (presumably themselves) to
 * allow merging. */

template<class Gene, template<class> class Derived>
class Visitor {

  using SP = std::shared_ptr<Gene>;

protected:

  virtual SP visit(const SP& self, const Derived<Gene>&) {
    return self;
  }

}; // class Visitor


/* Chain template dependency is used to generate the visit() methods one for
 * each of a list of derived classes. */

template<class Gene,
  template<class> class Head,
  template<class> class... Tail>
class Visitors :
  public Visitor<Gene, Head>,
  public Visitors<Gene, Tail...>
{

public:

  using Visitor<Gene, Head>::visit;
  using Visitors<Gene, Tail...>::visit;

}; // class Visitors

template<class Gene,
  template<class> class Last>
class Visitors<Gene, Last> :
  public Visitor<Gene, Last>
{

public:

  using Visitor<Gene, Last>::visit;

}; // class Visitors<Gene, Last>


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

} // namespace QGA
