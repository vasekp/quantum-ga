/* Enable or disable colours:
 * Colours::use = isatty(1); */

namespace Colours {

namespace {

  bool use{false};

  constexpr const char* BOLD = "\033[1m";
  constexpr const char* RED = "\033[1;31m";
  constexpr const char* GREEN = "\033[1;32m";
  constexpr const char* YELLOW = "\033[1;33m";
  constexpr const char* BLUE = "\033[1;34m";
  constexpr const char* RESET = "\033[0m";

} // anonymous inner namespace


template<typename... T>
class ColourPrinter {

public:

  ColourPrinter(const char* col_, const T&... ref_):
    col(col_), tuple(ref_...) { }

  friend std::ostream& operator<< (std::ostream& os, const ColourPrinter& cp) {
    if(use)
      os << cp.col;
    cp.write(os);
    if(use)
      os << Colours::RESET;
    return os;
  }

private:

  template<size_t Index = 0>
  void write(std::ostream& os) const {
    os << std::get<Index>(tuple);
    constexpr bool next = Index + 1 < sizeof...(T);
    if(next)
      // Conditioning on a known boolean is ugly but better still than
      // referring to an incomplete type if Index + 1 is too large
      write<next ? Index+1 : 0>(os);
  }

  const char* col;
  std::tuple<const T&...> tuple;

}; // class ColourPrinter<T...>


template<typename... T>
ColourPrinter<T...> bold(const T&... ref) {
  return {BOLD, ref...};
}

template<typename... T>
ColourPrinter<T...> red(const T&... ref) {
  return {RED, ref...};
}

template<typename... T>
ColourPrinter<T...> green(const T&... ref) {
  return {GREEN, ref...};
}

template<typename... T>
ColourPrinter<T...> yellow(const T&... ref) {
  return {YELLOW, ref...};
}

template<typename... T>
ColourPrinter<T...> blue(const T&... ref) {
  return {BLUE, ref...};
}

} // namespace Colours
