namespace QGA {

namespace internal {
  // Defined below
  template<typename...>
  class DomTuple;

  template<class...>
  class Counter;
}


template<class... Elements>
struct Fitness : public internal::DomTuple<Elements...> {

  using internal::DomTuple<Elements...>::DomTuple;

}; // struct Fitness<Elements...>

template<typename... Gates, typename Main, typename... Others>
struct Fitness<internal::Counter<Gates...>, Main, Others...> {

  internal::DomTuple<Main, Others...> tuple;
  internal::Counter<Gates...> counter;

  // Lexicographical ordering (first element of fitness most important)
  friend bool operator< (const Fitness& a, const Fitness& b) {
    return a.tuple < b.tuple
      || (a.tuple == b.tuple && a.counter < b.counter);
  }

  friend bool operator<< (const Fitness& a, const Fitness& b) {
    return (a.tuple <<= b.tuple) && (a.counter <<= b.counter) && !(a == b);
  }

  friend bool operator== (const Fitness& a, const Fitness& b) {
    return a.tuple == b.tuple && a.counter == b.counter;
  }

  friend double dist(const Fitness& a, const Fitness& b) {
    return dist(a.tuple, b.tuple) + dist(a.counter, b.counter);
  }

  operator Main() const {
    return static_cast<Main>(tuple);
  }

  friend std::ostream& operator<< (std::ostream& os, const Fitness& f) {
    return os << '{' << f.tuple << ',' << f.counter << '}';
  }

  friend std::istream& operator>> (std::istream& is, Fitness& f) {
    return is >> f.tuple >> f.counter;
  }

}; // struct Fitness<Counter, Elements...>


namespace internal {

/* Helper template class for DomTuple and Counter. Holds one element of a
 * specified type and defines operators for dominance comparison,
 * lexicoraphical comparison, and for output. Inherited in a zig-zag pattern,
 * for example:
 * Counter<A, B, C>
 *   DomComparator<int, Counter<B, C>>
 * Counter<B, C>
 *   DomCompataror<int, Counter<C>>
 * Counter<C> (specialization)
 *   DomComparator<int, void) (specialization)
 */

template<typename Element, class Next>
class DomComparator : protected Next {

public:
  template<typename... Args>
  DomComparator(Element elm_, Args... args): Next(args...), element(elm_) { }

  DomComparator() = default;

  friend bool operator<<= (const DomComparator& c1, const DomComparator& c2) {
    return c1.element <= c2.element && (c1.next() <<= c2.next());
  }

  friend bool operator== (const DomComparator& c1, const DomComparator& c2) {
    return c1.element == c2.element && c1.next() == c2.next();
  }

  friend bool operator< (const DomComparator& c1, const DomComparator& c2) {
    return c1.element < c2.element
      || (c1.element == c2.element && c1.next() < c2.next());
  }

  friend double dist(const DomComparator& c1, const DomComparator& c2) {
    return std::abs(c1.element - c2.element) + dist(c1.next(), c2.next());
  }

  operator const Element&() const {
    return element;
  }

  friend std::ostream& operator<< (std::ostream& os, const DomComparator& c) {
    return os << c.element << ',' << c.next();
  }

  friend std::istream& operator>> (std::istream& is, DomComparator& c) {
    is >> c.element;
    if(!is) {
      is.clear(is.rdstate() & ~std::ios::failbit);
      std::string str;
      is >> str;
      c.element = (Element)INFINITY;
    }
    return is >> c.next();
  }

protected:

  operator Element&() {
    return element;
  }

private:

  Next& next() {
    return static_cast<Next&>(*this);
  }

  const Next& next() const {
    return static_cast<const Next&>(*this);
  }

  Element element{};

}; // class DomComparator<Element, Next>

template<typename Element>
class DomComparator<Element, void> {

public:

  DomComparator(Element elm_): element(elm_) { }

  DomComparator() = default;

  friend bool operator<<= (const DomComparator& c1, const DomComparator& c2) {
    return c1.element <= c2.element;
  }

  friend bool operator== (const DomComparator& c1, const DomComparator& c2) {
    return c1.element == c2.element;
  }

  friend bool operator< (const DomComparator& c1, const DomComparator& c2) {
    return c1.element < c2.element;
  }

  friend double dist(const DomComparator& c1, const DomComparator& c2) {
    return std::abs(c1.element - c2.element);
  }

  operator const Element&() const {
    return element;
  }

  friend std::ostream& operator<< (std::ostream& os, const DomComparator& c) {
    return os << c.element;
  }

  friend std::istream& operator>> (std::istream& is, DomComparator& c) {
    is >> c.element;
    if(!is) {
      is.clear(is.rdstate() & ~std::ios::failbit);
      std::string str;
      is >> str;
      c.element = (Element)INFINITY;
    }
    return is;
  }

protected:

  operator Element&() {
    return element;
  }

private:

  Element element{};

}; // class DomComparator<Element, void> (sequence terminator)


/* Specialized as DomTuple<A, B, C, ...>, this template holds one content
 * element of each type given in the parameter pack. The types are expected to
 * be numeric and can repeat. It defines equality and dominance comparison
 * operators and a stream operator<< needed for a seamless inclusion in
 * Fitness. It is also implicitly convertible to A&, yielding a reference to
 * the first element. */

template<typename Head, typename... Tail>
class DomTuple<Head, Tail...> : public DomComparator<Head, DomTuple<Tail...>> {

public:

  using DomComparator<Head, DomTuple<Tail...>>::DomComparator;

  friend std::istream& operator>> (std::istream& is, DomTuple& c) {
    return is >> static_cast<DomComparator<Head, DomTuple<Tail...>>&>(c);
  }

}; // class DomTuple<Head, Tail...>

template<typename Last>
class DomTuple<Last> : public DomComparator<Last, void> {

public:

  using DomComparator<Last, void>::DomComparator;

  friend std::istream& operator>> (std::istream& is, DomTuple& c) {
    return is >> static_cast<DomComparator<Last, void>&>(c);
  }

}; // class DomTuple<Last> (sequence terminator)


/* Specialized as Counter<A, B, C, ...>, this template holds one counter per
 * each type given in the parameter pack. The types must be unique. It defines
 * member functions hit(const A*), hit(const B*), ...  which each bump their
 * respective counter.  Otherwise it behaves like a DomTuple<int, int, ...>. */

template<class Head, class... Tail>
class Counter<Head, Tail...> :
  public DomComparator<int, Counter<Tail...>>
{

public:

  void hit(const Head*) {
    ++(int&)(*this);
  }

  using Counter<Tail...>::hit;

  friend std::istream& operator>> (std::istream& is, Counter& c) {
    return is >> static_cast<DomComparator<int, Counter<Tail...>>&>(c);
  }

}; // class Counter<Head, Tail...>

template<class Last>
class Counter<Last> : public DomComparator<int, void> {

public:

  void hit(const Last*) {
    ++(int&)(*this);
  }

  friend std::istream& operator>> (std::istream& is, Counter& c) {
    return is >> static_cast<DomComparator<int, void>&>(c);
  }

}; // class Counter<Last> (sequence terminator)

} // namespace internal

} // namespace QGA
