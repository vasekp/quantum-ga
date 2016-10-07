namespace QGA {


template<class, template<class> class, template<class> class...>
class Visitors;


/* The base class for all genes. Defines methods derived classes have to
 * implement, and provides default (no-op) definition for some of them. */
template<class Gene,
  template<class> class... Derived>
class GeneBase : public Visitors<Gene, Derived...> {

  using SP = std::shared_ptr<Gene>;

public:

  // apply this gene to a state vector
  virtual Wrapper::State apply(const Wrapper::State&) const = 0;

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

  static unsigned ctrlBitString(unsigned rand, unsigned target) {
    unsigned ctrl = 0;
    double c = (double)rand / std::numeric_limits<unsigned>::max();
    /* Convert an unsigned between 0 and UINT_MAX to a bit string where the
     * probability of 1 in each position is given by Config::pControl. A value
     * less than 0.5 means that plain NOTs and C-NOTs will be generated more
     * often than CC-NOTs and higher. */
    for(unsigned i = 0; i < Config::nBit - 1; i++) {
      ctrl <<= 1;
      if(c < Config::pControl) {
        ctrl |= 1;
        c /= Config::pControl;
      } else {
        c = (c - Config::pControl)/(1 - Config::pControl);
      }
    }
    /* At this point ctrl has nBit - 1 bits. We use this to guarantee that
     * 1 << target is left unoccupied. */
    return
      ((ctrl >> target) << (target + 1))  // shift bits left of tgt to the left
        |
      (ctrl & ((1 << target) - 1));  // keep bits right of tgt
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


/* This is the default Gene class. It can be used when the operations above
 * are sufficient and the gene polymorphism reduces to choosing a subset from
 * the existing Gene classes. A different class, however, must be defined when
 * additional functionality is required in the base (see wrapper-search for an
 * example). */

template<template<class> class... Derived>
class Gene: public GeneBase<Gene<Derived...>, Derived...> {

public:

  static std::shared_ptr<Gene> getNew();

}; // class Gene

} // namespace QGA
