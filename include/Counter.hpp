namespace QGA {

namespace internal {

/* This template accepts a series of typenames and holds a corresponding
 * number of zero-initialized counters. It defines functions hit(Type_k*) for
 * each type which bump the respective counter. */

template<class GateBase, class Head, class... Tail>
class Counter: public Counter<GateBase, Tail...> {

  unsigned cnt = 0;

public:

  Counter() = default;

  void hit(const typename Head::template Template<GateBase>*) {
    cnt++;
  }

  using Counter<GateBase, Tail...>::hit;

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

  const Counter<GateBase, Tail...>& next() const {
    return static_cast<const Counter<GateBase, Tail...>&>(*this);
  }

}; // class Counter<GateBase, Head, Tail...>


template<class GateBase, class Last>
class Counter<GateBase, Last> {

  unsigned cnt = 0;

public:

  Counter() = default;

  void hit(const typename Last::template Template<GateBase>*) {
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

}; // class Counter<GateBase, Last>

} // namespace internal

} // namespace QGA
