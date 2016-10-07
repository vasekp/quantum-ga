namespace QGA {


template<class, template<class> class, template<class> class...>
class Visitors;


template<class Gene,
  template<class> class... Derived>
class GeneBase : public Visitors<Gene, Derived...> {

public:

  virtual Wrapper::State apply(const Wrapper::State&) const = 0;

  virtual unsigned complexity() const = 0;

  virtual bool invert() {
    return false;
  }

  virtual bool mutate() {
    return false;
  }

  virtual bool simplify() {
    return false;
  }

  /* http://www.oodesign.com/visitor-pattern.html */
  bool merge(const Gene* next) {
    return next->invite(static_cast<Gene*>(this));
  }

  friend std::ostream& operator<< (std::ostream& os, const GeneBase& g) {
    return g.write(os);
  }

protected:

  virtual bool invite(Gene*) const = 0;

  virtual std::ostream& write(std::ostream&) const = 0;

}; // virtual class GeneBase<Gene, Derived...>


template<class Gene, template<class> class Derived>
class Visitor {

protected:

  virtual bool visit(const Derived<Gene>&) {
    return false;
  }

};


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

};

template<class Gene,
  template<class> class Last>
class Visitors<Gene, Last> :
  public Visitor<Gene, Last>
{

public:

  using Visitor<Gene, Last>::visit;

};


template<template<class> class... Derived>
class Gene: public GeneBase<Gene<Derived...>, Derived...> {

public:

  static Gene* getNew();

};

} // namespace QGA
