namespace QGA {

namespace internal {

/* Called as Chooser<Base, A, B, C, ...>::getNew(i) returns:
 * i = 0: A<Base>::getNew()
 * i = 1: B<Base>::getNew()
 * i = 2: C<Base>::getNew()
 * etc. */

template<class GeneBase,
  template<class> class Head,
  template<class> class... Tail>
class Chooser {

public:

  static std::shared_ptr<GeneBase> getNew(unsigned index) {
    if(index == 0)
      return Head<GeneBase>::getNew();
    else
      return Chooser<GeneBase, Tail...>::getNew(index - 1);
  }

}; // class Chooser<GeneBase, Head, Tail...>

template<class GeneBase,
  template<class> class Last>
class Chooser<GeneBase, Last> {

public:

  static std::shared_ptr<GeneBase> getNew(unsigned index) {
    if(index == 0)
      return Last<GeneBase>::getNew();
    else
      throw std::logic_error("Index too large in Chooser!");
  }

}; // class Chooser<GeneBase, Last>

} // namespace internal


/* Given QGA::GeneBase or its subclass of the same template parameters,
 * along with a selection of gene templates, construct a gene ready for use
 * with a Candidate. Implements a static getNew() function which randomly
 * picks from the given Genes. */

template<template<class, template<class> class...> class GeneBase,
  template<class> class... Genes>
class CustomGene : public GeneBase<CustomGene<GeneBase, Genes...>, Genes...> {

public:

  static std::shared_ptr<CustomGene> getNew() {
    std::uniform_int_distribution<> dist(0, sizeof...(Genes) - 1);
    return internal::Chooser<CustomGene, Genes...>::getNew(dist(gen::rng));
  }

}; // class CustomGene<GeneBase, Genes...>


/* This is the default case with GeneBase = QGA::GeneBase for brevity. */

template<template<class> class... Genes>
using Gene = CustomGene<GeneBase, Genes...>;

} // namespace QGA
