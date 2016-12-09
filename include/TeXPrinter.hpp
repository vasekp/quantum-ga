class TeXPrinter : public QGA::CircuitPrinter {

public:

  TeXPrinter(unsigned nBit_): nBit(nBit_), lines{} {
    for(unsigned line = 0; line < nBit; line++)
      lines.push_back({});
  }

  /* Prints a single gate. */
  void addGate(std::string name, unsigned line) override {
    _addControlled(name, line, {});
  }

  /* Prints a single qubit controlled gate. */
  void addControlledGate(std::string name,
      unsigned line, std::vector<unsigned> controls) override {
    _addControlled(name, line, controls);
  }

  /* Prints a swap gate. */
  void addSwapGate(unsigned line1, unsigned line2) override {
    _addSwap(line1, line2);
  }

  /* Prints a gate spanning all lines of output. */
  void addBarrierGate(std::string name) override {
    _addBroadGate(0, nBit - 1, {name});
  }

  std::ostream& print(std::ostream& os) const override {
    alignAll();
    for(auto& line : lines) {
      for(auto& element : line)
        os << element;
      os << "& \\qw \\\\\n";
    }
    return os;
  }

private:

  void _addControlled(const std::string& name,
      unsigned line, const std::vector<unsigned>& controls) {
    unsigned first, last;

    /* Find the total span of the gates being printed */
    if(controls.size() > 0) {
      auto minmax = std::minmax_element(controls.begin(), controls.end());
      first = *minmax.first;
      last = *minmax.second;
      if(line < first)
        first = line;
      if(line > last)
        last = line;
    } else
      first = last = line;

    align(first, last);

    for(unsigned i = first, prev = first; i <= last; i++) {
      if(i == line) {
        std::ostringstream os{};
        os << "& \\gate{" << name << "} \\qwx[-" << i - prev << "] ";
        lines[i].push_back(os.str());
        prev = i;
      } else {
        auto it = std::find(controls.begin(), controls.end(), i);
        if(it != controls.end()) { // i is one of the control qubits
          std::ostringstream os{};
          os << "& \\ctrl{-" << i - prev << "} ";
          lines[i].push_back(os.str());
          prev = i;
        }
      }
    }

    align(first, last);
  }

  void _addSwap(unsigned line1, unsigned line2) {
    if(line1 > line2)
      std::swap(line1, line2);
    align(line1, line2);
    {
      std::ostringstream os{};
      os << "& \\qswap \\qwx[" << line2 - line1 << "] ";
      lines[line1].push_back(os.str());
    }
    lines[line2].push_back("& \\qswap ");
    align(line1, line2);
  }

  void _addBroadGate(unsigned lineFrom, unsigned lineTo, std::string name) {
    alignAll();
    {
      std::ostringstream os{};
      os << "& \\multigate{" << lineTo - lineFrom << "}{" << name << "} ";
      lines[lineFrom].push_back(os.str());
    }
    for(unsigned i = lineFrom + 1; i <= lineTo; i++)
      lines[i].push_back("& \\ghost{" + name + "} ");
  }

  void align(unsigned lineFrom, unsigned lineTo) const {
    /* Find the longest line wihin the span so far */
    size_t max = 0;
    for(unsigned i = lineFrom; i <= lineTo; i++)
      if(lines[i].size() > max)
        max = lines[i].size();

    /* Pad the lines with a qubit wire */
    for(unsigned i = lineFrom; i <= lineTo; i++)
      for(size_t end = lines[i].size(); end < max; end++)
        lines[i].push_back("& \\qw ");
  }

  void alignAll() const {
    align(0, nBit - 1);
  }

  unsigned nBit;
  // alignment may be needed from print which takes a const reference
  // NB that aligning does not technically change the contents
  mutable std::vector<std::vector<std::string>> lines;

}; // class TeXPrinter
