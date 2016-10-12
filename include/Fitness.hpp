namespace QGA {

struct Fitness {

  double error;
  size_t length;
  unsigned cplx;

  friend std::ostream& operator<< (std::ostream& os, const Fitness& f) {
    return os << '{'
       << f.error << ','
       << f.length << ','
       << f.cplx << '}';
  }

  friend bool operator< (const Fitness& a, const Fitness& b) {
    return a.error < b.error;
  }

  friend bool operator<< (const Fitness& a, const Fitness& b) {
    return a.error <= b.error
        && a.length <= b.length
        && a.cplx <= b.cplx
        && !(a == b);
  }

  friend bool operator== (const Fitness& a, const Fitness& b) {
    return a.error == b.error
        && a.length == b.length
        && a.cplx == b.cplx;
  }

}; // struct Fitness

} // namespace QGA
