#ifndef REGEX_HPP
#define REGEX_HPP

#include <string>
#include <memory>
#include <stddef.h>

/* The purpose of these classes is to allow compiling the std::regex template
 * sources separately and then use the functionality via an opaque library, as
 * recompiling them every time with -O3 is very time-consuming. */

namespace regex {

  class matches {

    class matches_impl;

  public:

    matches();
    ~matches();
    bool matched(size_t index);
    std::string match(size_t index);
    matches_impl& impl();

  private:

    std::unique_ptr<matches_impl> pImpl;

  }; // class matches


  class regex {

    class regex_impl;

  public:

    regex(std::string expr);
    ~regex();
    bool match(std::string searched, matches& ms);

  private:

    std::unique_ptr<regex_impl> pImpl;

  }; // class regex

} // namespace regex

#endif // !defined REGEX_HPP
