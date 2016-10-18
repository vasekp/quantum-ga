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

} // namespace QGA
