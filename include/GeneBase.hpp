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

  /* The following functions pass a std::shared_ptr (SP) pointing to this
   * along with this. If a gene allows a given operation, it can use this
   * parameter to actually rewrite the shared pointer by a new
   * std::make_shared<>(...). If it does not, it can return with a no-op.
   * The internals of SP guarantee that in either case all reference counts
   * are updated properly and no memory is leaked.
   *
   * IMPORTANT: reassigning self from within the functions can result in a
   * deletion of this. Keep a local copy of the shared pointer on stack if
   * that is anything than the very last command. */

  virtual void invert(SP& /*self*/) { }

  virtual void mutate(SP& /*self*/) { }

  virtual void simplify(SP& /*self*/) { }

  /* Merge needs to be implemented using double dispatch, because only
   * genes of the same class can typically be merged (albeit this is not a
   * restriction) but all we know at compile time is a pointer to the base
   * class.
   *
   * Let first is a shared pointer to Derived1 and second is a shared pointer
   * to Derived2. Both classes override their invite() and visit() methods.
   * The call pattern goes as following:
   *   first->merge(first, second)
   *   -> second->invite(first, second) [vtable lookup in Derived2]
   *   -> first->visit(first, second, Derived2&) [vtable lookup in Derived1]
   * Now the override Derived1::visit() is called in an overload with its
   * third parameter being a (const) Derived2&. This allows to react properly
   * to any possible combination of the two derived classes.
   *
   * See: http://www.oodesign.com/visitor-pattern.html */

  bool merge(SP& first, SP& second) {
    return second->invite(first, second);
  }

  friend std::ostream& operator<< (std::ostream& os, const GeneBase& g) {
    return g.write(os);
  }

protected:

  /* Every derived class must implement this function with exactly the
   * following definition:
   *
   *   bool invite(SP& first, SP& second) const override {
   *     return first->visit(first, second, *this);
   *   }
   *
   * This can't be done here because *this only refers to the derived class
   * itself in its own context. We need to call a particular visitor for a
   * specific class so it can't be a const GeneBase&. See the description of
   * merge() above.
   *
   * The derived class (gene template) needs to define this function even if
   * it does not allow merging with any other genes. */

  virtual bool invite(SP& first, SP& second) const = 0;

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


/* The purpose of this helper class is to inject a virtual method for calling
 * a particular gene class, for the purposes of the merge() call pattern. Note
 * that for the visitor design pattern, we need one method like this for each
 * possible visitee, but we don't know what the derived classes are yet. (This
 * is supplied in Gene.hpp.)
 *
 * Subclasses can override visit(SP&, SP&, const X&) for some particular
 * values of X (presumably themselves) to allow merging with instances of
 * class X.
 *
 * The default implementation returns false, indicating that no merge
 * happened. A return value of true means that one of the pair of genes has
 * been consumed. Note that a particular implementation can also return false
 * if it actually modifies first and second but does not combine them into a
 * single gene. */

template<class Gene, template<class> class Derived>
class Visitor {

  using SP = std::shared_ptr<Gene>;

protected:

  virtual bool visit(SP& /*first*/, SP& /*second*/, const Derived<Gene>&) {
    return false;
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
