namespace QGA {

class CircuitPrinter {

public:

  /* Prints a single gate. */
  virtual void addGate(const std::string& name, unsigned line) = 0;

  /* Prints a single qubit controlled gate. */
  virtual void addControlledGate(const std::string& name,
      unsigned line, const std::vector<unsigned>& controls) = 0;

  /* Prints a swap gate. */
  virtual void addSwapGate(unsigned line1, unsigned line2) = 0;

  /* Prints a gate spanning all lines of output. */
  virtual void addBarrierGate(const std::string& name) = 0;

  friend std::ostream& operator<< (std::ostream& os, const CircuitPrinter& p) {
    return p.print(os);
  }

  virtual ~CircuitPrinter() { };

private:

  virtual std::ostream& print(std::ostream& os) const = 0;

}; // class CircuitPrinter

} // namespace QGA
