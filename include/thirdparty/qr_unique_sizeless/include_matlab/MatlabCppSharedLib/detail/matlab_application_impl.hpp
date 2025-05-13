/* Copyright 2017 The MathWorks, Inc. */

#ifndef MATLAB_APPLICATION_IMPL_HPP
#define MATLAB_APPLICATION_IMPL_HPP

#include "../cppsharedlib_api.hpp"
#include "../matlab_application.hpp"

namespace matlab {
namespace cpplib {

inline MATLABApplication::MATLABApplication(
    const MATLABApplicationMode mode,
    const std::vector<std::u16string> &options)
    : _mode(mode), _options(options)
{
}

inline MATLABApplication::~MATLABApplication()
{
  terminated = true;
  runtime_terminate_session();
}

} // namespace cpplib
} // namespace matlab

#endif /* MATLAB_APPLICATION_IMPL_HPP */
