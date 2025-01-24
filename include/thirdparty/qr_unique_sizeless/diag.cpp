//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// diag.cpp
//
// Code generation for function 'diag'
//

// Include files
#include "diag.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
namespace coder {
void diag(const array<double, 2U> &v, array<double, 1U> &d)
{
  if ((v.size(0) == 1) && (v.size(1) == 1)) {
    d.set_size(1);
    d[0] = v[0];
  } else {
    int u0;
    int u1;
    u0 = v.size(0);
    u1 = v.size(1);
    if (u0 <= u1) {
      u1 = u0;
    }
    if (v.size(1) > 0) {
      u0 = u1;
    } else {
      u0 = 0;
    }
    d.set_size(u0);
    for (u1 = 0; u1 < u0; u1++) {
      d[u1] = v[u1 + v.size(0) * u1];
    }
  }
}

} // namespace coder

// End of code generation (diag.cpp)
