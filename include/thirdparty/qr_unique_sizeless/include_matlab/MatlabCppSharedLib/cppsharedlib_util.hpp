/* Copyright 2017 The MathWorks, Inc. */

#ifndef CPPSHAREDLIB_UTIL_HPP
#define CPPSHAREDLIB_UTIL_HPP

#include "cppsharedlib_api.hpp"

#include <algorithm>
#include <cstring>
#include <mutex>
#include <ratio>

namespace matlab {
namespace cpplib {

class MATLABLibrary;

enum class MATLABApplicationMode
{
  OUT_OF_PROCESS = 0,
  IN_PROCESS = 1
};
} // namespace cpplib
} // namespace matlab

#endif // CPPSHAREDLIB_UTIL_HPP