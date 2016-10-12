namespace QGA {

/* Forward declarations */
namespace Backend {
  class State;
}

namespace internal {
  template<class, template<class> class, template<class> class...>
  class Visitors;
}
/* End forward declarations */


/* The base class for all genes. Defines methods derived classes have to
 * implement, and provides default (no-op) definition for some of them. */
template<class RealBase,
  template<class> class... Gates>
class GateBase : public internal::Visitors<RealBase, Gates...> {

public:

  using Pointer = std::shared_ptr<RealBase>;

  // apply this gene to a state vector
  virtual Backend::State applyTo(const Backend::State&) const = 0;

  // return an arbitrary notion of complexity of this operation (accumulative)
  virtual unsigned complexity() const = 0;

  // return whether this gate has degenerated to the identity (e.g., by means
  // of simplification or merge)
  virtual bool isTrivial() {
    return false;
  }

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

  virtual void invert(Pointer& /*self*/) { }

  virtual void mutate(Pointer& /*self*/) { }

  virtual void simplify(Pointer& /*self*/) { }

  /* Merge needs to be implemented using double dispatch, because only
   * genes of the same class can typically be merged (albeit this is not a
   * restriction) but all we know at compile time is a pointer to the base
   * class.
   *
   * Let first is a shared pointer to Derived1 and second is a shared pointer
   * to Derived2. Both classes override their invite() and 3-argument merge()
   * methods. The call pattern goes as following:
   *   first->merge(first, second)
   *   -> second->invite(first, second) [vtable lookup in Derived2]
   *   -> first->merge(first, second, Derived2&) [vtable lookup in Derived1]
   * Now the override Derived1::merge() is called in an overload with its
   * third parameter being a (const) Derived2&. This allows to react properly
   * to any possible combination of the two derived classes.
   *
   * See: http://www.oodesign.com/visitor-pattern.html */

  bool merge(Pointer& first, Pointer& second) {
    if(first->isTrivial()) {
      // op1 = identity: replace by second and consume
      first = second;
      return true;
    } else if(second->isTrivial()) {
      // op2 = identity: consume
      return true;
    } else
      return second->invite(first, second);
  }

  using internal::Visitors<RealBase, Gates...>::merge;

  friend std::ostream& operator<< (std::ostream& os, const GateBase& g) {
    return g.write(os);
  }

protected:

  /* Every derived class must implement this function with exactly the
   * following definition:
   *
   *   bool invite(Pointer& first, Pointer& second) const override {
   *     return first->merge(first, second, *this);
   *   }
   *
   * This can't be done here because *this only refers to the derived class
   * itself in its own context. We need to call a particular visitor for a
   * specific class so it can't be a const GateBase&. See the description of
   * merge() above.
   *
   * The derived class (gene template) needs to define this function even if
   * it does not allow merging with any other genes. */

  virtual bool invite(Pointer& first, Pointer& second) const = 0;

  virtual std::ostream& write(std::ostream&) const = 0;

}; // virtual class GateBase<RealBase, Gates...>


namespace internal {

/* The purpose of this helper class is to inject a virtual method for calling
 * a particular gene class, for the purposes of the merge() call pattern. Note
 * that for the visitor design pattern, we need one method like this for each
 * possible visitee, but we don't know what the derived classes are yet. (This
 * is supplied in Gene.hpp.)
 *
 * Subclasses can override merge(Pointer&, Pointer&, const X&) for some
 * particular values of X (presumably themselves) to allow merging with
 * instances of class X.
 *
 * The default implementation returns false, indicating that no merge
 * happened. A return value of true means that one of the pair of genes has
 * been consumed. Note that a particular implementation can also return false
 * if it actually modifies first and second but does not combine them into a
 * single gene. */

template<class GateBase, template<class> class Gate>
class Visitor {

  using SP = std::shared_ptr<GateBase>;

protected:

  virtual bool merge(SP& /*first*/, SP& /*second*/, const Gate<GateBase>&) {
    return false;
  }

}; // class Visitor


/* Chain template dependency is used to generate the merge() methods one for
 * each of a list of derived classes. */

template<class GateBase,
  template<class> class Head,
  template<class> class... Tail>
class Visitors :
  public Visitor<GateBase, Head>,
  public Visitors<GateBase, Tail...>
{

public:

  using Visitor<GateBase, Head>::merge;
  using Visitors<GateBase, Tail...>::merge;

}; // class Visitors

template<class GateBase,
  template<class> class Last>
class Visitors<GateBase, Last> :
  public Visitor<GateBase, Last>
{

public:

  using Visitor<GateBase, Last>::merge;

}; // class Visitors<GateBase, Last>

} // namespace internal

} // namespace QGA
