namespace QGA {

template<class Counter>
struct Fitness {

  double error;
  unsigned cplx;
  Counter cc;

  friend std::ostream& operator<< (std::ostream& os, const Fitness& f) {
    return os << '{'
       << f.error << ','
       << f.cplx << ','
       << f.cc << '}';
  }

  friend bool operator< (const Fitness& a, const Fitness& b) {
    return a.error < b.error;
  }

  friend bool operator<< (const Fitness& a, const Fitness& b) {
    return a.error <= b.error
        && a.cplx <= b.cplx
        && (a.cc <<= b.cc)
        && !(a == b);
  }

  friend bool operator== (const Fitness& a, const Fitness& b) {
    return a.error == b.error
        && a.cplx == b.cplx
        && (a.cc == b.cc);
  }

}; // struct Fitness


namespace internal {

/* Specialized as Counter<A, B, C, ...>, this template holds one counter per
 * each type given in the parameter pack. It defines member functions
 * hit(const A*), hit(const B*), ... which each bump their respective counter.
 * There are also dominance comparison operators and a stream operator<<
 * needed for a seamless inclusion in Fitness. */

template<class Head, class... Tail>
class Counter: public Counter<Tail...> {

  unsigned cnt = 0;

public:

  Counter() = default;

  void hit(const Head*) {
    cnt++;
  }

  using Counter<Tail...>::hit;

  friend std::ostream& operator<< (std::ostream& os, const Counter& c) {
    return os << c.cnt << ',' << c.next();
  }

  friend bool operator<<= (const Counter& c1, const Counter& c2) {
    return c1.cnt <= c2.cnt && (c1.next() <<= c2.next());
  }

  friend bool operator== (const Counter& c1, const Counter& c2) {
    return c1.cnt == c2.cnt && c1.next() == c2.next();
  }

  friend bool operator<< (const Counter& c1, const Counter& c2) {
    return c1 <<= c2 && !(c1 == c2);
  }

private:

  const Counter<Tail...>& next() const {
    return static_cast<const Counter<Tail...>&>(*this);
  }

}; // class Counter<Head, Tail...>

template<class Last>
class Counter<Last> {

  unsigned cnt = 0;

public:

  Counter() = default;

  void hit(const Last*) {
    cnt++;
  }

  friend std::ostream& operator<< (std::ostream& os, const Counter& c) {
    return os << c.cnt;
  }

  friend bool operator<<= (const Counter& c1, const Counter& c2) {
    return c1.cnt <= c2.cnt;
  }

  friend bool operator== (const Counter& c1, const Counter& c2) {
    return c1.cnt == c2.cnt;
  }

  friend bool operator<< (const Counter& c1, const Counter& c2) {
    return c1 < c2;
  }

}; // class Counter<Last>

} // namespace internal

} // namespace QGA
