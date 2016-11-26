#ifndef QGA_HPP
#define QGA_HPP

#include "QGA_commons.hpp"
#include "genetic.hpp"
#include "regex.hpp"

#include "QGA_bits/Backend.hpp"
#include "QGA_bits/CircuitPrinter.hpp"
#include "QGA_bits/Fitness.hpp"
#include "QGA_bits/GateBase.hpp" // uses Fitness.hpp and CircuitPrinter.hpp
#include "QGA_bits/Gene.hpp"     // uses GateBase.hpp
#include "QGA_bits/Tools.hpp"
#include "QGA_bits/Gates.hpp"    // uses Tools.hpp and GateBase.hpp
#include "QGA_bits/CandidateCounter.hpp"
#include "QGA_bits/CandidateBase.hpp"
#include "QGA_bits/CandidateFactory.hpp"

#endif // !defined QGA_HPP
