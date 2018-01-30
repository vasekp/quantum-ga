namespace QGA {

namespace internal {
  // Defined below
  template<class, class, class...>
  class CastInject;
}


/* The base class for all gates. Defines methods derived classes have to
 * implement, and provides default (no-op) definition for some of them. */

template<class ContextParm, class... Gates>
class GateBase :
  protected internal::CastInject<GateBase<ContextParm, Gates...>, Gates...>
{

protected:

  using Context = ContextParm;

public:

  // Both required in QGA::Gene
  using Pointer = std::shared_ptr<const GateBase>;
  using Counter = internal::Counter<
    typename Gates::template Template<GateBase>...
  >;

  // apply this gate to a state vector
  virtual Backend::State applyTo(const Backend::State&,
      const Context* = nullptr) const = 0;

  // return the number of control qubits of this gate
  virtual unsigned controls() const {
    return 0;
  }

  /* The following functions pass a std::shared_ptr (SP) pointing to this
   * along with this. If a gate allows a given operation, it should return a
   * new SP created using std::make_shared<OwnClass>. If it does not, it
   * should return the parameter self so it can be reused just with its
   * reference count upped instead of creating a new copy. The default
   * implementation does precisely that. */

  virtual Pointer getAnother() const = 0;

  virtual Pointer invert(const Pointer& self) const {
    return self;
  }

  virtual Pointer mutate(const Pointer& self) const {
    return self;
  }

  virtual Pointer simplify(const Pointer& self) const {
    return self;
  }

  virtual Pointer swapQubits(const Pointer&, unsigned, unsigned) const = 0;

  /* The cast() function is a light-weight version of dynamic_cast: knowing
   * all the derived subclasses *a priori*, we can define a virtual function
   * for each that returns a nullptr if the type does not match and returns
   * this (casted to the particular derived type) if it does. The former is
   * the default behaviour inherited from CastInject but the latter must be
   * explicitly implemented in each subclass like this:
   *
   *   const Derived* cast(const Derived*) const override {
   *     return this;
   *   }
   */
  using internal::CastInject<GateBase, Gates...>::cast;

  virtual bool sameType(const GateBase&) const {
    return false;
  }

  // returns {} (NULL shared pointer) if merge can't be done
  virtual Pointer merge(const GateBase&) const = 0;

  // return whether this gate has degenerated to the identity (e.g., by means
  // of simplification or merge)
  virtual bool isTrivial() const {
    return false;
  }

  // Called from CandidateBase
  virtual void hit(Counter& c) const = 0;

  friend std::ostream& operator<< (std::ostream& os, const GateBase& g) {
    return g.write(os);
  }

  // print this gate in circuit representation
  virtual void printOn(CircuitPrinter& p) const = 0;

  virtual ~GateBase() { }

private:

  virtual std::ostream& write(std::ostream&) const = 0;

  /* Finally, the derived classes need to provide the following static method:
   *
   *   static Pointer read(const std::string& s);
   *
   * This can not be enforced or default-implemented by an abstract base
   * class. A no-op can return Pointer{}. */

}; // virtual class GateBase<Context, Gates...>


namespace internal {

/* The purpose of this helper class is to inject a virtual method
 * cast(const Derived*). Note that we need one method like this for each
 * possible derived class. The idea is that all functions return a nullptr
 * except the one that corresponds to the actual derived class. That must be
 * overridden in that class's definition to return this. */

template<class GateBase, class Gate>
class SingleCastInject {

  using Pointer = std::shared_ptr<const GateBase>;
  using GateResolved = typename Gate::template Template<GateBase>;

public:

  virtual const GateResolved* cast(const GateResolved*) const {
    return nullptr;
  }

  virtual ~SingleCastInject() { }

}; // class SingleCastInject


/* Chain template dependency is used to generate the cast() methods one for
 * each of a parameter pack of derived classes (gate templates).  The list of
 * possible gate types to be considered is specific to each problem and is
 * provided through a specialization of QGA::Gene. */

template<class GateBase, class Head, class... Tail>
class CastInject :
  SingleCastInject<GateBase, Head>,
  CastInject<GateBase, Tail...>
{

public:

  using SingleCastInject<GateBase, Head>::cast;
  using CastInject<GateBase, Tail...>::cast;

}; // class CastInject

template<class GateBase, class Last>
class CastInject<GateBase, Last> : SingleCastInject<GateBase, Last> {

public:

  using SingleCastInject<GateBase, Last>::cast;

}; // class CastInject<GateBase, Last>

} // namespace internal

} // namespace QGA
