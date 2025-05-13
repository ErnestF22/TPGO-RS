/* Copyright 2017-2020 The MathWorks, Inc. */

#ifndef ENGINE_UTIL_HPP
#define ENGINE_UTIL_HPP

#include "MatlabExecutionInterface/util.hpp"
#include "cpp_engine_api.hpp"

#include <algorithm>
#include <cstring>
#include <mutex>
#include <ratio>

namespace matlab {

namespace engine {

typedef matlab::execution::StreamBuffer StreamBuffer;
typedef std::u16string String;
class MATLABEngine;

enum class WorkspaceType
{
  BASE = 0,
  GLOBAL = 1
};
} // namespace engine
} // namespace matlab

#endif // ENGINE_UTIL_HPP