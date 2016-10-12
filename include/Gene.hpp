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

  static std::shared_ptr<GateBase> getNew(unsigned index) {
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

  static std::shared_ptr<GateBase> getNew(unsigned index) {
    if(index == 0)
      return Last<GateBase>::getNew();
    else
      throw std::logic_error("Index too large in Chooser!");
  }

}; // class Chooser<GateBase, Last>

} // namespace internal


/* Given QGA::GateBase or its subclass of the same template parameters,
 * along with a selection of gene templates, construct a gene ready for use
 * with a Candidate. Implements a static getNew() function which randomly
 * picks from the given Genes. */

template<template<class, template<class> class...> class GateBase,
  template<class> class... Genes>
class CustomGene : public GateBase<CustomGene<GateBase, Genes...>, Genes...> {

public:

  static std::shared_ptr<CustomGene> getNew() {
    std::uniform_int_distribution<> dist(0, sizeof...(Genes) - 1);
    return internal::Chooser<CustomGene, Genes...>::getNew(dist(gen::rng));
  }

}; // class CustomGene<GateBase, Genes...>


/* This is the default case with GateBase = QGA::GateBase for brevity. */

template<template<class> class... Genes>
using Gene = CustomGene<GateBase, Genes...>;

} // namespace QGA
