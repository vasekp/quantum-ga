namespace QGA {

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

public:

  static typename GateBase::Pointer getNew(unsigned index) {
    if(index == 0)
      return Head<GateBase>::getNew();
    else
      return Chooser<GateBase, Tail...>::getNew(index - 1);
  }

}; // class Chooser<GateBase, Head, Tail...>

template<class GateBase,
  template<class> class Last>
class Chooser<GateBase, Last> {

public:

  static typename GateBase::Pointer getNew(unsigned index) {
    if(index == 0)
      return Last<GateBase>::getNew();
    else
      throw std::logic_error("Index too large in Chooser!");
  }

}; // class Chooser<GateBase, Last>

} // namespace internal


template<template<class, template<class> class...> class GateBase,
  template<class> class... Genes>
class Gate : public GateBase<Gate<GateBase, Genes...>, Genes...> { };


/* Given QGA::GateBase or its subclass of the same template parameters,
 * along with a selection of gene templates, construct a gene ready for use
 * with a Candidate.
 *
 * This implements a static getNew() function which randomly picks from the
 * given Genes, along with several shortcuts to functions of the GateBase.
 * Other functions are redirected using a ->. */

template<template<class, template<class> class...> class GateBase,
  template<class> class... Genes>
class CustomGene : Gate<GateBase, Genes...>::Pointer {

  using CGate = Gate<GateBase, Genes...>;
  using Pointer = typename CGate::Pointer;

public:

  CustomGene(const Pointer& ptr): Pointer(ptr) { }

  CustomGene(Pointer&& ptr): Pointer(std::move(ptr)) { }

  static Pointer getNew() {
    std::uniform_int_distribution<> dist(0, sizeof...(Genes) - 1);
    return internal::Chooser<CGate, Genes...>::getNew(dist(gen::rng));
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

}; // class CustomGene<GateBase, Genes...>


/* This is the default case with GateBase = QGA::GateBase for brevity. */

template<template<class> class... Genes>
using Gene = CustomGene<GateBase, Genes...>;

} // namespace QGA
