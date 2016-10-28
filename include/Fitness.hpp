namespace QGA {

template<class Counter>
struct Fitness {

  double error;
  Counter cc;
  unsigned controls;

  friend std::ostream& operator<< (std::ostream& os, const Fitness& f) {
    return os << '{'
       << f.error << ','
       << f.cc << ','
       << f.controls << '}';
  }

  friend bool operator< (const Fitness& a, const Fitness& b) {
    return a.error < b.error;
  }

  friend bool operator<< (const Fitness& a, const Fitness& b) {
    return a.error <= b.error
        && (a.cc <<= b.cc)
        && a.controls <= b.controls
        && !(a == b);
  }

  friend bool operator== (const Fitness& a, const Fitness& b) {
    return a.error == b.error
        && (a.cc == b.cc)
        && a.controls == b.controls;
  }

}; // struct Fitness


namespace internal {

/* Specialized as DomTuple<A, B, C, ...>, this template holds one content
 * element of each type given in the parameter pack. The types are expected to
 * be numeric and can repeat. It defines equality and dominance comparison
 * operators and a stream operator<< needed for a seamless inclusion in
 * Fitness. */

template<typename... Types>
class DomTuple;

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

template<class...>
class Counter;

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
