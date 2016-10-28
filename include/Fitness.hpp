namespace QGA {

namespace internal {
  // Defined below
  template<typename...>
  class DomTuple;

  template<class...>
  class Counter;
}


template<class... Elements>
struct Fitness : internal::DomTuple<Elements...> {

  using internal::DomTuple<Elements...>::DomTuple;

}; // struct Fitness<Elements...>

template<typename... Gates, typename Main, typename... Others>
struct Fitness<internal::Counter<Gates...>, Main, Others...> {

  internal::DomTuple<Main, Others...> tuple;
  internal::Counter<Gates...> counter;

  friend std::ostream& operator<< (std::ostream& os, const Fitness& f) {
    return os << '{' << f.tuple << ',' << f.counter << '}';
  }

  // Comparison of the first element of tuple
  friend bool operator< (const Fitness& a, const Fitness& b) {
    return static_cast<Main>(a.tuple) < static_cast<Main>(b.tuple);
  }

  friend bool operator<< (const Fitness& a, const Fitness& b) {
    return (a.tuple <<= b.tuple) && (a.counter <<= b.counter) && !(a == b);
  }

  friend bool operator== (const Fitness& a, const Fitness& b) {
    return a.tuple == b.tuple && a.counter == b.counter;
  }

  operator Main() const {
    return static_cast<Main>(tuple);
  }

}; // struct Fitness<Counter, Elements...>


namespace internal {

/* Specialized as DomTuple<A, B, C, ...>, this template holds one content
 * element of each type given in the parameter pack. The types are expected to
 * be numeric and can repeat. It defines equality and dominance comparison
 * operators and a stream operator<< needed for a seamless inclusion in
 * Fitness. */

template<typename Head, typename... Tail>
class DomTuple<Head, Tail...> : virtual DomTuple<Tail...> {

public:

  DomTuple(Head head, Tail... tail) :
    DomTuple<Tail...>(tail...), element(head) { }

  friend std::ostream& operator<< (std::ostream& os, const DomTuple& c) {
    os << c.element;
    if(sizeof...(Tail))
      os << ',' << c.next();
    return os;
  }

  friend bool operator<<= (const DomTuple& c1, const DomTuple& c2) {
    return c1.element <= c2.element && (c1.next() <<= c2.next());
  }

  friend bool operator== (const DomTuple& c1, const DomTuple& c2) {
    return c1.element == c2.element && c1.next() == c2.next();
  }

  operator Head() const {
    return element;
  }

  template<typename New>
  using Prepend = DomTuple<New, Head, Tail...>;

protected:

  // Default initialization is only intended for derived classes (Counter)
  DomTuple(): element() { }

  Head element;

private:

  const DomTuple<Tail...>& next() const {
    return static_cast<const DomTuple<Tail...>&>(*this);
  }

}; // class DomTuple<Head, Tail...>

template<>
class DomTuple<> {

public:

  // no-op
  friend std::ostream& operator<< (std::ostream& os, const DomTuple&) {
    return os;
  }

  // trivial
  friend bool operator<<= (const DomTuple&, const DomTuple&) {
    return true;
  }

  // trivial
  friend bool operator== (const DomTuple&, const DomTuple&) {
    return true;
  }

  template<typename New>
  using Prepend = DomTuple<New>;

}; // class DomTuple<> (trivial)


/* Specialized as Counter<A, B, C, ...>, this template holds one counter per
 * each type given in the parameter pack. The types must be unique and
 * non-void. It defines member functions hit(const A*), hit(const B*), ...
 * which each bump their respective counter.  Otherwise it behaves like a
 * DomTuple<unsigned, unsigned, ...>. */

template<class Head, class... Tail>
class Counter<Head, Tail...> :
  Counter<Tail...>,
  virtual public Counter<Tail...>::Tuple::template Prepend<unsigned>
{

public:

  using Tuple = typename Counter<Tail...>::Tuple::template Prepend<unsigned>;

  void hit(const Head*) {
    Tuple::element++;
  }

  using Counter<Tail...>::hit;

}; // class Counter<Head, Tail...>

template<>
class Counter<> : virtual DomTuple<> {

public:

  using Tuple = DomTuple<>;

  // unimplemented: not callable
  void hit(void*) = delete;

}; // class Counter<> (trivial)

} // namespace internal

} // namespace QGA
