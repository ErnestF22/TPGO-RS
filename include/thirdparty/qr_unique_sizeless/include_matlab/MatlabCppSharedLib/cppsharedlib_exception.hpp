/* Copyright 2017 The MathWorks, Inc. */

#ifndef CPPSHAREDLIB_EXCEPTION_HPP
#define CPPSHAREDLIB_EXCEPTION_HPP

#include <MatlabExecutionInterface/detail/exception_impl.hpp>
#include <MatlabExecutionInterface/exception.hpp>

namespace matlab {

namespace cpplib {

class CppSharedLibException final : public matlab::execution::Exception {
public:
  CppSharedLibException();

  /**
   * Constructor that accepts an error message
   */
  CppSharedLibException(const std::string &msg);
};

} // namespace cpplib
} // namespace matlab

#endif /* CPPSHAREDLIB_EXCEPTION_HPP */