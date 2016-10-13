namespace QGA {

namespace internal {
  // Defined below
  template<class, template<class> class, template<class> class...>
  class Chooser;
}


/* Given QGA::GateBase or its subclass of the same template parameters,
 * along with a selection of gene templates, construct a gene ready for use
 * with a CandidateBase.
 *
 * This holds a shared pointer to GateBase (as typedef'd in GateBase.hpp) and
 * as such can refer to different gates polymorphically.
 *
 * Implements a static getNew() function which randomly picks from the given
 * Gates, along with several shortcuts to functions of the GateBase.  Other
 * functions (like complexity() or functions added by extensions to the
 * original QGA::GateBase) are redirected using a ->. */

template<template<class, template<class> class...> class GateBase,
  template<class> class... Gates>
class CustomGene : internal::Gate<GateBase, Gates...>::Pointer {

  using CGate = internal::Gate<GateBase, Gates...>;
  using Pointer = typename CGate::Pointer;

public:

  CustomGene(const Pointer& ptr): Pointer(ptr) { }

  CustomGene(Pointer&& ptr): Pointer(std::move(ptr)) { }

  static Pointer getNew() {
    return internal::Chooser<CGate, Gates...>::getNew();
  }

  const CGate* operator->() const {
    return Pointer::operator->();
  }

  friend std::ostream& operator<<(std::ostream& os, const CustomGene& g) {
    return os << *g;
  }

  /* These helper functions simplify the call pattern of the pointer-passing
   * virtual functions of GateBase. */
  void invert() {
    pointer()->invert(pointer());
  }

  void mutate() {
    pointer()->mutate(pointer());
  }

  void simplify() {
    pointer()->simplify(pointer());
  }

  bool merge(CustomGene& other) {
    return pointer()->merge(pointer(), other);
  }

  /* Two genes are equal iff they point to the same object. This is used in
   * CandidateFactory to check if a genotype has changed at all during a
   * mutation (if not, the parent is returned). */
  bool operator==(const CustomGene& other) const {
    return pointer() == other.pointer();
  }

private:

  Pointer& pointer() {
    return static_cast<Pointer&>(*this);
  }

  const Pointer& pointer() const {
    return static_cast<const Pointer&>(*this);
  }

}; // class CustomGene<GateBase, Gates...>


/* Unless the application intends to extend QGA::GateBase we use it as the
 * default parameter. It can then simply refer to QGA::Gene<Genes...>. If an
 * extended functionality is requested of the GateBase this becomes
 * QGA::CustomGene<NewBase, Genes...>. */

template<template<class> class... Gates>
using Gene = CustomGene<GateBase, Gates...>;


namespace internal {

/* Called as Chooser<Base, A, B, C, ...>::getNew(i) returns:
 * i = 0: A<Base>::getNew()
 * i = 1: B<Base>::getNew()
 * i = 2: C<Base>::getNew()
 * etc. */

template<class GateBase,
  template<class> class Head,
  template<class> class... Tail>
class Chooser {

  template<class, template<class> class, template<class> class...>
  friend class Chooser;

  static typename GateBase::Pointer getNew(unsigned index) {
    if(index == 0)
      return Head<GateBase>::getNew();
    else
      return Chooser<GateBase, Tail...>::getNew(index - 1);
  }

public:

  static typename GateBase::Pointer getNew() {
    // upper bound inclusive: no need to add 1 for Head
    std::uniform_int_distribution<> dist(0, sizeof...(Tail));
    return getNew(dist(gen::rng));
  }

}; // class Chooser<GateBase, Head, Tail...>

template<class GateBase,
  template<class> class Last>
class Chooser<GateBase, Last> {

public:

#ifdef DEBUG
  static typename GateBase::Pointer getNew(unsigned index = 0) {
    if(index == 0)
      return Last<GateBase>::getNew();
    else
      throw std::logic_error("Index too large in Chooser!");
  }
#else
  static typename GateBase::Pointer getNew(unsigned = 0) {
    return Last<GateBase>::getNew();
  }
#endif

}; // class Chooser<GateBase, Last>

} // namespace internal

} // namespace QGA
