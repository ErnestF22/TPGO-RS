/* Copyright 2017 The MathWorks, Inc. */

#ifndef ENGINE_EXCEPTION_IMPL_HPP
#define ENGINE_EXCEPTION_IMPL_HPP

#include "../engine_exception.hpp"
#include <future>
#include <memory>
#include <streambuf>
#include <string>
#include <vector>

#if defined(_WIN32)
#define NOEXCEPT throw()
#else
#define NOEXCEPT noexcept
#endif

namespace matlab {
namespace engine {

inline EngineException::EngineException()
{
}

inline EngineException::EngineException(const std::string &msg) : Exception(msg)
{
}

} // namespace engine
} // namespace matlab

#endif /* ENGINE_EXCEPTION_IMPL_HPP */