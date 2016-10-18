namespace QGA {

namespace internal {
  // Defined below
  template<class, class, class...>
  class Visitors;
}


/* The base class for all genes. Defines methods derived classes have to
 * implement, and provides default (no-op) definition for some of them. */

template<class ContextParm, class... Gates>
class GateBase :
  public internal::Visitors<GateBase<ContextParm, Gates...>, Gates...>
{

public:

  using Pointer = std::shared_ptr<const GateBase>;
  using Counter = internal::Counter<
    typename Gates::template Template<GateBase>...>;
  using Context = ContextParm;

  // apply this gene to a state vector
  virtual Backend::State applyTo(const Backend::State&, const Context*) const
    = 0;

  // return an arbitrary notion of complexity of this operation (accumulative)
  virtual unsigned complexity() const = 0;

  // return whether this gate has degenerated to the identity (e.g., by means
  // of simplification or merge)
  virtual bool isTrivial() const {
    return false;
  }

  /* The following functions pass a std::shared_ptr (SP) pointing to this
   * along with this. If a gene allows a given operation, it should return a
   * new SP created using std::make_shared<OwnClass>. If it does not, it
   * should return the parameter self so it can be reused just with its
   * reference count upped instead of creating a new copy. The default
   * implementation does precisely that. */

  virtual Pointer invert(const Pointer& self) const {
    return self;
  }

  virtual Pointer mutate(const Pointer& self) const {
    return self;
  }

  virtual Pointer simplify(const Pointer& self) const {
    return self;
  }

  /* Merge needs to be implemented using double dispatch, because only
   * genes of the same class can typically be merged (albeit this is not a
   * restriction) but all we know at compile time is a pointer to the base
   * class.
   *
   * Let first is a shared pointer to Derived1 and second is a shared pointer
   * to Derived2. Both classes override their invite() and 1-argument merge()
   * methods. The call pattern goes as following:
   *   first->merge(first, second)
   *   -> second->invite(first)         [vtable lookup in Derived2]
   *   -> first->merge(const Derived2&) [vtable lookup in Derived1]
   * Now the override Derived1::merge() is called in an overload with its
   * third parameter being a const Derived2&. This allows to react properly
   * to any possible combination of the two derived classes.
   *
   * See: http://www.oodesign.com/visitor-pattern.html */

  Pointer merge(const Pointer& first, const Pointer& second) const {
    if(first->isTrivial()) {
      // op1 = identity: consume and return other
      return second;
    } else if(second->isTrivial()) {
      // op2 = identity: consume and return first
      return first;
    } else
      return second->invite(first);
  }

  using internal::Visitors<GateBase, Gates...>::merge;

  /* This function allows us to count the number of each gate type using
   * QGA::internal::Counter. However each gate must call it itself so that the
   * Counter can recognize the type the hit request came from. Therefore each
   * gate must implement this verbatim:
   *
   *   void hit(typename GateBase::Counter& c) const {
   *     c.hit(this);
   *   }
   */

  virtual void hit(Counter& c) const = 0;

  friend std::ostream& operator<< (std::ostream& os, const GateBase& g) {
    return g.write(os);
  }

  virtual ~GateBase() { }

protected:

  /* Every derived class must implement this function with exactly the
   * following definition:
   *
   *   Pointer invite(const Pointer& first) const override {
   *     return first->merge(*this);
   *   }
   *
   * This can't be done here because *this only refers to the derived class
   * itself in its own context. We need to call a particular visitor for a
   * specific class so it can be recognized as a (const) GateBase&. See the
   * description of merge() above.
   *
   * The derived class (gene template) needs to define this function even if
   * it does not allow merging with any other genes. */

  virtual Pointer invite(const Pointer& other) const = 0;

  virtual std::ostream& write(std::ostream&) const = 0;

  /* Finally, the derived classes need to provide the following static method:
   *
   *   static Pointer read(const std::string& s);
   *
   * This can not be enforced or default-implemented by an abstract base
   * class. A no-op can return Pointer{}. */

}; // virtual class GateBase<Context, Gates...>


namespace internal {

/* The purpose of this helper class is to inject a virtual method for calling
 * a particular gene class, for the purposes of the merge() call pattern. Note
 * that for the visitor design pattern, we need one method like this for each
 * possible visitee, but we don't know what the derived classes are yet. (This
 * is supplied in Gene.hpp.)
 *
 * Subclasses can override merge(const X&) for some particular values of X
 * (presumably themselves) to allow merging with instances of class X.
 *
 * The default implementation returns an unassigned std::shared_pointer,
 * indicating that no merge happened. Any other return value means that the
 * pair of genes has been consumed and should be replaced by the return value
 * (which can also be one of them). */

template<class GateBase, class Gate>
class Visitor {

  using Pointer = std::shared_ptr<const GateBase>;
  using GateResolved = typename Gate::template Template<GateBase>;

protected:

  virtual Pointer merge(const GateResolved&) const {
    return {};
  }

  virtual ~Visitor() { }

}; // class Visitor


/* Chain template dependency is used to generate the merge() methods one for
 * each of a list of derived classes. */

template<class GateBase, class Head, class... Tail>
class Visitors :
  public Visitor<GateBase, Head>,
  public Visitors<GateBase, Tail...>
{

public:

  using Visitor<GateBase, Head>::merge;
  using Visitors<GateBase, Tail...>::merge;

}; // class Visitors

template<class GateBase, class Last>
class Visitors<GateBase, Last> :
  public Visitor<GateBase, Last>
{

public:

  using Visitor<GateBase, Last>::merge;

}; // class Visitors<GateBase, Last>

} // namespace internal

} // namespace QGA
