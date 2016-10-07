namespace QGA {


template<class, template<class> class, template<class> class...>
class Visitors;


template<class Gene,
  template<class> class... Derived>
class GeneBase : public Visitors<Gene, Derived...> {

  using SP = std::shared_ptr<Gene>;

public:

  virtual Wrapper::State apply(const Wrapper::State&) const = 0;

  virtual unsigned complexity() const = 0;

  virtual SP invert(const SP& self) {
    return self;
  }

  virtual SP mutate(const SP& self) {
    return self;
  }

  virtual SP simplify(const SP& self) {
    return self;
  }

  /* http://www.oodesign.com/visitor-pattern.html */
  SP merge(const SP& self, const SP& other) {
    return other.get()->invite(self);
  }

  friend std::ostream& operator<< (std::ostream& os, const GeneBase& g) {
    return g.write(os);
  }

protected:

  virtual SP invite(const SP&) const = 0;

  virtual std::ostream& write(std::ostream&) const = 0;

}; // virtual class GeneBase<Gene, Derived...>


template<class Gene, template<class> class Derived>
class Visitor {

  using SP = std::shared_ptr<Gene>;

protected:

  virtual SP visit(const SP& self, const Derived<Gene>&) {
    return self;
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

  static std::shared_ptr<Gene> getNew();

};

} // namespace QGA
