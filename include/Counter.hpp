namespace QGA {

namespace internal {

/* This template accepts a series of typenames and holds a corresponding
 * number of zero-initialized counters. It defines functions hit(Type_k*) for
 * each type which bump the respective counter. */

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
