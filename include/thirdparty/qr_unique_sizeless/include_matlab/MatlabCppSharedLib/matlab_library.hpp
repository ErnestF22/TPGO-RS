/* Copyright 2017-2022 The MathWorks, Inc. */

#ifndef MATLAB_LIBRARY_HPP
#define MATLAB_LIBRARY_HPP

#include "matlab_application.hpp"
#include <MatlabExecutionInterface/execution_interface.hpp>
#include <complex>
#include <future>
#include <memory>
#include <streambuf>
#include <vector>

namespace matlab {

namespace cpplib {

using namespace matlab::execution;
class MATLABLibrary : public matlab::execution::ExecutionInterface {
public:
  /**
   * Destructor
   *
   * @throw none
   */
  ~MATLABLibrary();

  /**
   * wait for all figures to be closed
   **/
  void waitForFiguresToClose();

private:
  friend FutureResult<std::unique_ptr<MATLABLibrary>>
  initMATLABLibraryAsync(std::shared_ptr<MATLABApplication> application,
                         const std::u16string &ctffilename,
                         const std::u16string &session_key);
  /**
   * Constructor
   *
   * @param handle - The internal implementation
   *
   * @throw none
   */
  MATLABLibrary(std::shared_ptr<MATLABApplication> application,
                uint64_t handle);

  std::shared_ptr<MATLABApplication> appPtr_;
};

} // namespace cpplib
} // namespace matlab

#endif // MATLAB_LIBRARY_HPP
