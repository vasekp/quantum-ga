namespace QGA {

/* Specialized as Fitness<A, B, C, ...>, this template holds one content
 * element of each type given in the parameter pack. The types are expected to
 * be numeric and can repeat. The class defines equality and dominance
 * comparison operators, a distance measure, and stream input and output
 * operators. */

template<typename...>
class Fitness;

template<typename Head, typename... Tail>
class Fitness<Head, Tail...> : protected Fitness<Tail...> {

  using Next = Fitness<Tail...>;
  using Element = Head;
  
  template<typename...>
  friend class Fitness;

public:
  template<typename... Args>
  Fitness(Element elm_, Args... args): Next(args...), element(elm_) { }

  Fitness() = default;

  friend bool operator<<= (const Fitness& c1, const Fitness& c2) {
    return c1.element <= c2.element && (c1.next() <<= c2.next());
  }

  friend bool operator<< (const Fitness& c1, const Fitness& c2) {
    return (c1 <<= c2) && !(c1 == c2);
  }

  friend bool operator== (const Fitness& c1, const Fitness& c2) {
    return c1.element == c2.element && c1.next() == c2.next();
  }

  friend bool operator< (const Fitness& c1, const Fitness& c2) {
    return c1.element < c2.element
      || (c1.element == c2.element && c1.next() < c2.next());
  }

  friend double dist(const Fitness& c1, const Fitness& c2) {
    return (c1.element > c2.element
      ? c1.element - c2.element
      : c2.element - c1.element) + dist(c1.next(), c2.next());
  }

  const Element& head() const {
    return element;
  }

  friend std::ostream& operator<< (std::ostream& os, const Fitness& c) {
    os << '{';
    c.format(os);
    return os << '}';
  }

  friend std::istream& operator>> (std::istream& is, Fitness& c) {
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

  void format(std::ostream& os) const {
    os << element;
    if(sizeof...(Tail) > 0) {
      os << ',';
      next().format(os);
    }
  }

private:

  Next& next() {
    return static_cast<Next&>(*this);
  }

  const Next& next() const {
    return static_cast<const Next&>(*this);
  }

  Element element{};

}; // class Fitness<Head, Tail>


// Stub class for sequence termination
template<>
class Fitness<> {

public:

  Fitness() = default;

  friend bool operator<<= (const Fitness&, const Fitness&) {
    return true;
  }

  friend bool operator<< (const Fitness&, const Fitness&) {
    return true;
  }

  friend bool operator== (const Fitness&, const Fitness&) {
    return true;
  }

  friend bool operator< (const Fitness&, const Fitness&) {
    return false;
  }

  friend double dist(const Fitness&, const Fitness&) {
    return 0;
  }

  void format(std::ostream&) const { }

  friend std::istream& operator>> (std::istream& is, Fitness&) {
    return is;
  }

}; // class Fitness<> (sequence terminator)

} // namespace QGA
