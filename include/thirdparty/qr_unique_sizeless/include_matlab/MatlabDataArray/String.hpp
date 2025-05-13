/* Copyright 2016-2017 The MathWorks, Inc. */

#ifndef MATLAB_DATA_STRING_HPP
#define MATLAB_DATA_STRING_HPP

#include "Optional.hpp"
#include <string>

namespace matlab {
namespace data {
using String = std::basic_string<char16_t>;

using MATLABString = optional<String>;
} // namespace data
} // namespace matlab

#endif
