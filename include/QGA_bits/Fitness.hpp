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

  friend std::ostream& operator<< (std::ostream& os, const Fitness& f) {
    return os << '{' << f.tuple << ',' << f.counter << '}';
  }

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

}; // struct Fitness<Counter, Elements...>


namespace internal {

/* Helper template class for DomTuple and Counter. Holds one element of a
 * specified type and defines operators for dominance comparison,
 * lexicoraphical comparison, and for output. Inherited in a zig-zag pattern,
 * for example:
 * Counter<A, B, C>
 *   DomComparator<unsigned, Counter<B, C>>
 * Counter<B, C>
 *   DomCompataror<unsigned, Counter<C>>
 * Counter<C> (specialization)
 *   DomComparator<unsigned, void) (specialization)
 */

template<typename Element, class Next>
class DomComparator : protected Next {

public:
  template<typename... Args>
  DomComparator(Element elm_, Args... args): Next(args...), element(elm_) { }

  friend std::ostream& operator<< (std::ostream& os, const DomComparator& c) {
    return os << c.element << ',' << c.next();
  }

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

protected:

  operator Element&() {
    return element;
  }

  DomComparator() = default;

private:

  const Next next() const {
    return static_cast<const Next&>(*this);
  }

  Element element{};

}; // class DomComparator<Element, Next>

template<typename Element>
class DomComparator<Element, void> {

public:

  DomComparator(Element elm_): element(elm_) { }

  friend std::ostream& operator<< (std::ostream& os, const DomComparator& c) {
    return os << c.element;
  }

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

protected:

  operator Element&() {
    return element;
  }

  DomComparator() = default;

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

  DomTuple(Head head, Tail... tail) :
    DomComparator<Head, DomTuple<Tail...>>(head, tail...) { }

}; // class DomTuple<Head, Tail...>

template<typename Last>
class DomTuple<Last> : public DomComparator<Last, void> {

public:

  DomTuple(Last last) :
    DomComparator<Last, void>(last) { }

}; // class DomTuple<Last> (sequence terminator)


/* Specialized as Counter<A, B, C, ...>, this template holds one counter per
 * each type given in the parameter pack. The types must be unique. It defines
 * member functions hit(const A*), hit(const B*), ...  which each bump their
 * respective counter.  Otherwise it behaves like a DomTuple<unsigned,
 * unsigned, ...>. */

template<class Head, class... Tail>
class Counter<Head, Tail...> :
  public DomComparator<unsigned, Counter<Tail...>>
{

public:

  void hit(const Head*) {
    ++(unsigned&)(*this);
  }

  using Counter<Tail...>::hit;

}; // class Counter<Head, Tail...>

template<class Last>
class Counter<Last> : public DomComparator<unsigned, void> {

public:

  void hit(const Last*) {
    ++(unsigned&)(*this);
  }

}; // class Counter<Last> (sequence terminator)

} // namespace internal

} // namespace QGA
