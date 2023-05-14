#pragma once

#ifndef PSEUDO_FLOAT_IOSTREAM_CPP_H
#define PSEUDO_FLOAT_IOSTREAM_CPP_H

// Use a separate file for this, as don't want to always have iostream included to use this class

#include "PseudoFloat.h"
#include <iostream>

// Create a better version of this later if needed. The hack is to just to convert it to a double
inline std::ostream &operator<<(std::ostream &os, PseudoFloat const &x) { 
    return os << (double)x;
}
#endif
