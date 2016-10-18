namespace QGA {

namespace internal {
  // Defined below
  template<class, class...>
  class Chooser;

  template<class, class...>
  class Reader;
}


/* The Gene type ready for use with a CandidateBase.
 *
 * This holds a shared pointer to GateBase (as typedef'd in GateBase.hpp) and
 * as such can refer to different gates polymorphically.
 *
 * Implements a static getNew() function which randomly picks from the given
 * Gates, along with several shortcuts to functions of the GateBase.  Other
 * functions like applyTo() or hit() are redirected using a ->. */

template<class Context, class... Gates>
class CustomGene : GateBase<Context, Gates...>::Pointer {

  using CGate = GateBase<Context, Gates...>;
  using Pointer = typename CGate::Pointer;

public:

  using Counter = typename CGate::Counter;

  CustomGene() = default; // Needed in CandidateBase::read()

  CustomGene(const Pointer& ptr): Pointer(ptr) { }

  CustomGene(Pointer&& ptr): Pointer(std::move(ptr)) { }

  static Pointer getNew() {
    return internal::Chooser<
      typename Gates::template Template<CGate>...
    >::getNew();
  }

  const CGate* operator->() const {
    return Pointer::operator->();
  }

  friend std::ostream& operator<<(std::ostream& os, const CustomGene& g) {
    return os << *g;
  }

  friend std::istream& operator>>(std::istream& is, CustomGene& g) {
    std::string gene{};
    if(!(is >> gene))
      return is;
    Pointer ptr = internal::Reader<
      typename Gates::template Template<CGate>...
    >::read(gene);
    if(!ptr) {
      is.setstate(std::ios::failbit);
      return is;
    }
    g = {ptr};
    return is;
  }

  /* These helper functions simplify the call pattern of the pointer-passing
   * virtual functions of GateBase. */
  void invert() {
    pointer() = pointer()->invert(pointer());
  }

  void mutate() {
    pointer() = pointer()->mutate(pointer());
  }

  void simplify() {
    pointer() = pointer()->simplify(pointer());
  }

  bool merge(CustomGene& other) {
    Pointer result = pointer()->merge(pointer(), other);
    if(result) {
      pointer() = result;
      return true;
    } else
      return false;
  }

  /* Two genes are equal iff they point to the same object. This is used in
   * CandidateFactory to check if a genotype has changed at all during a
   * mutation (if not, the parent is returned). */
  bool operator==(const CustomGene& other) const {
    return pointer() == other.pointer();
  }

private:

  Pointer& pointer() {
    return static_cast<Pointer&>(*this);
  }

  const Pointer& pointer() const {
    return static_cast<const Pointer&>(*this);
  }

}; // class CustomGene<Context, Gates...>


/* The default value for Context (additional data passed to GateBase::applyTo)
 * is void. Unfortunately there's no way of providing default values for
 * template parameters before the parameter pack without invalidating the use
 * of the pack so we use a template class alias instead. */

template<class... Gates>
using Gene = CustomGene<void, Gates...>;


namespace internal {

/* Called as Chooser<A, B, C, ...>::getNew(i), returns:
 * i = 0: A::getNew()
 * i = 1: B::getNew()
 * i = 2: C::getNew()
 * etc. */

template<class Head, class... Tail>
class Chooser {

  template<class, class...>
  friend class Chooser;

  static auto getNew(unsigned index) -> decltype(Head::getNew()) {
    if(index == 0)
      return Head::getNew();
    else
      return Chooser<Tail...>::getNew(index - 1);
  }

public:

  static auto getNew() -> decltype(Head::getNew()) {
    // upper bound inclusive: no need to add 1 for Head
    std::uniform_int_distribution<> dist(0, sizeof...(Tail));
    return getNew(dist(gen::rng));
  }

}; // class Chooser<Head, Tail...>

template<class Last>
class Chooser<Last> {

public:

#ifdef DEBUG
  static auto getNew(unsigned index = 0) -> decltype(Last::getNew()) {
    if(index == 0)
      return Last::getNew();
    else
      throw std::logic_error("Index too large in Chooser!");
  }
#else
  static auto getNew(unsigned = 0) -> decltype(Last::getNew()) {
    return Last::getNew();
  }
#endif

}; // class Chooser<Last>


/* Called as Reader<A, B, C, ...>::read(input), this helper tries calling
 * A::read(input), B::read(input), C::read(input)... until the first returns a
 * value convertible to boolean true. */

template<class Head, class... Tail>
class Reader {

public:

  static auto read(std::string input) -> decltype(Head::read(input)) {
    if(auto ret = Head::read(input))
      return ret;
    else
      return Reader<Tail...>::read(input);
  }

}; // class Reader<Head, Tail...>

template<class Last>
class Reader<Last> {

public:

  static auto read(std::string input) -> decltype(Last::read(input)) {
    return Last::read(input);
  }

}; // class Reader<Last>

} // namespace internal

} // namespace QGA
