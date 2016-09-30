#include <unistd.h>

/* Enable or disable colours:
 * Colours::use = isatty(1); */

namespace Colours {

  bool use{false};

  const char* bold() { return use ? "\033[1m" : ""; }
  const char* red() { return use ? "\033[1;31m" : ""; }
  const char* green() { return use ? "\033[1;32m" : ""; }
  const char* yellow() { return use ? "\033[1;33m" : ""; }
  const char* blue() { return use ? "\033[1;34m" : ""; }
  const char* reset() { return use ? "\033[0m" : ""; }

} // namespace Colours
