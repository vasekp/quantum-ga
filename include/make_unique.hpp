#ifndef MAKE_UNIQUE_HPP
#define MAKE_UNIQUE_HPP

/* Not defined by the standard until C++14 */
template<class T, typename... Args>
static std::unique_ptr<T> make_unique(Args... args) {
  return std::unique_ptr<T>{new T(std::forward<Args>(args)...)};
}

#endif // !defined MAKE_UNIQUE_HPP
