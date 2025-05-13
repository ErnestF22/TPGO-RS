//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// qr_unique.cpp
//
// Code generation for function 'qr_unique'
//

// Include files
#include "qr_unique.h"
#include "coder_array.h"
#include "diag.h"
#include "qr.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <emmintrin.h>

// Function Declarations
static void binary_expand_op(coder::array<double, 2U> &in1,
                             const coder::array<double, 2U> &in2,
                             const coder::array<double, 1U> &in3);

static void binary_expand_op_1(coder::array<double, 2U> &in1,
                               const coder::array<double, 2U> &in2,
                               const coder::array<double, 1U> &in3);

// Function Definitions
static void binary_expand_op(coder::array<double, 2U> &in1,
                             const coder::array<double, 2U> &in2,
                             const coder::array<double, 1U> &in3)
{
  coder::array<double, 2U> b_in2;
  int b_loop_ub;
  int in3_idx_0;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  in3_idx_0 = in3.size(0);
  if (in3_idx_0 == 1) {
    loop_ub = in2.size(0);
  } else {
    loop_ub = in3_idx_0;
  }
  b_loop_ub = in2.size(1);
  b_in2.set_size(loop_ub, b_loop_ub);
  stride_0_0 = (in2.size(0) != 1);
  stride_1_0 = (in3_idx_0 != 1);
  for (int i{0}; i < b_loop_ub; i++) {
    for (int i1{0}; i1 < loop_ub; i1++) {
      b_in2[i1 + b_in2.size(0) * i] =
          in2[i1 * stride_0_0 + in2.size(0) * i] * in3[i1 * stride_1_0];
    }
  }
  in3_idx_0 = in1.size(0);
  stride_1_0 = in1.size(1);
  for (int i{0}; i < stride_1_0; i++) {
    for (int i1{0}; i1 < in3_idx_0; i1++) {
      in1[i1 + in1.size(0) * i] = b_in2[i1 + in3_idx_0 * i];
    }
  }
}

static void binary_expand_op_1(coder::array<double, 2U> &in1,
                               const coder::array<double, 2U> &in2,
                               const coder::array<double, 1U> &in3)
{
  coder::array<double, 2U> b_in2;
  int aux_0_1;
  int aux_1_1;
  int b_loop_ub;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  loop_ub = in2.size(0);
  if (in3.size(0) == 1) {
    b_loop_ub = in2.size(1);
  } else {
    b_loop_ub = in3.size(0);
  }
  b_in2.set_size(loop_ub, b_loop_ub);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_1 = (in3.size(0) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (int i{0}; i < b_loop_ub; i++) {
    int scalarLB;
    int vectorUB;
    scalarLB = (loop_ub / 2) << 1;
    vectorUB = scalarLB - 2;
    for (int i1{0}; i1 <= vectorUB; i1 += 2) {
      __m128d r;
      r = _mm_loadu_pd(&in2[i1 + in2.size(0) * aux_0_1]);
      _mm_storeu_pd(&b_in2[i1 + b_in2.size(0) * i],
                    _mm_mul_pd(r, _mm_set1_pd(in3[aux_1_1])));
    }
    for (int i1{scalarLB}; i1 < loop_ub; i1++) {
      b_in2[i1 + b_in2.size(0) * i] =
          in2[i1 + in2.size(0) * aux_0_1] * in3[aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  loop_ub = in1.size(0);
  b_loop_ub = in1.size(1);
  for (int i{0}; i < b_loop_ub; i++) {
    for (int i1{0}; i1 < loop_ub; i1++) {
      in1[i1 + in1.size(0) * i] = b_in2[i1 + loop_ub * i];
    }
  }
}

void qr_unique(const coder::array<double, 2U> &A, coder::array<double, 2U> &Q,
               coder::array<double, 2U> &R)
{
  __m128d b_r;
  coder::array<double, 2U> b_A;
  coder::array<double, 2U> q;
  coder::array<double, 2U> r;
  coder::array<double, 1U> s;
  int b_loop_ub;
  int k;
  int loop_ub;
  int nx;
  int scalarLB;
  int vectorUB;
  //  Thin QR factorization ensuring diagonal of R is real, positive if
  //  possible.
  //
  //  function [Q, R] = qr_unique(A)
  //
  //  If A is a matrix, then Q, R are matrices such that A = QR where Q'*Q = I
  //  and R is upper triangular. If A is real, then so are Q and R.
  //  This is a thin QR factorization in the sense that if A has more rows than
  //  columns, than Q has the same size as A.
  //
  //  If A has full column rank, then R has positive reals on its diagonal.
  //  Otherwise, it may have zeros on its diagonal.
  //
  //  This is equivalent to a call to Matlab's qr(A, 0), up to possible
  //  sign/phase changes of the columns of Q and of the rows of R to ensure
  //  the stated properties of the diagonal of R. If A has full column rank,
  //  this decomposition is unique, hence the name of the function.
  //
  //  If A is a 3D array, then Q, R are also 3D arrays and this function is
  //  applied to each slice separately.
  //
  //  See also: qr randrot randunitary
  //  This file is part of Manopt: www.manopt.org.
  //  Original author: Nicolas Boumal, June 18, 2019.
  //  Contributors:
  //  Change log:
  //    Sep. 24, 2023 (NB):
  //        Edited out bsxfun() for improved speed.
  if (A.size(0) >= A.size(1)) {
    //  A (or its slices) has more rows than columns
    Q.set_size(A.size(0), A.size(1));
    nx = A.size(0) * A.size(1);
    for (int i{0}; i < nx; i++) {
      Q[i] = 0.0;
    }
    R.set_size(A.size(1), A.size(1));
    nx = A.size(1) * A.size(1);
    for (int i{0}; i < nx; i++) {
      R[i] = 0.0;
    }
  } else {
    Q.set_size(A.size(0), A.size(0));
    nx = A.size(0) * A.size(0);
    for (int i{0}; i < nx; i++) {
      Q[i] = 0.0;
    }
    R.set_size(A.size(0), A.size(1));
    nx = A.size(0) * A.size(1);
    for (int i{0}; i < nx; i++) {
      R[i] = 0.0;
    }
  }
  loop_ub = A.size(0);
  b_loop_ub = A.size(1);
  b_A.set_size(A.size(0), A.size(1));
  for (int i{0}; i < b_loop_ub; i++) {
    for (int i1{0}; i1 < loop_ub; i1++) {
      b_A[i1 + b_A.size(0) * i] = A[i1 + A.size(0) * i];
    }
  }
  coder::qr(b_A, q, r);
  //  In the real case, s holds the signs of the diagonal entries of R.
  //  In the complex case, s holds the unit-modulus phases of these
  //  entries. In both cases, d = diag(s) is a unitary matrix, and
  //  its inverse is d* = diag(conj(s)).
  coder::diag(r, s);
  nx = s.size(0);
  for (k = 0; k < nx; k++) {
    if (std::isnan(s[k])) {
      s[k] = rtNaN;
    } else if (s[k] < 0.0) {
      s[k] = -1.0;
    } else {
      s[k] = (s[k] > 0.0);
    }
  }
  //  Since a = qr (with 'a' the slice of A currently being processed),
  //  it is also true that a = (qd)(d*r). By construction, qd still has
  //  orthonormal columns, and d*r has positive real entries on its
  //  diagonal, /unless/ s contains zeros. The latter can only occur if
  //  slice a does not have full column rank, so that the decomposition
  //  is not unique: we make an arbitrary choice in that scenario.
  //  While exact zeros are unlikely, they may occur if, for example,
  //  the slice a contains repeated columns, or columns that are equal
  //  to zero. If an entry should be mathematically zero but is only
  //  close to zero numerically, then it is attributed an arbitrary
  //  sign dictated by the numerical noise: this is also fine.
  nx = s.size(0);
  for (k = 0; k < nx; k++) {
    if (s[k] == 0.0) {
      s[k] = 1.0;
    }
  }
  if (s.size(0) == q.size(1)) {
    nx = q.size(0);
    k = q.size(1);
    b_A.set_size(q.size(0), q.size(1));
    for (int i{0}; i < k; i++) {
      scalarLB = (q.size(0) / 2) << 1;
      vectorUB = scalarLB - 2;
      for (int i1{0}; i1 <= vectorUB; i1 += 2) {
        b_r = _mm_loadu_pd(&q[i1 + q.size(0) * i]);
        _mm_storeu_pd(&b_A[i1 + b_A.size(0) * i],
                      _mm_mul_pd(b_r, _mm_set1_pd(s[i])));
      }
      for (int i1{scalarLB}; i1 < nx; i1++) {
        b_A[i1 + b_A.size(0) * i] = q[i1 + q.size(0) * i] * s[i];
      }
    }
    nx = A.size(0);
    k = Q.size(1);
    for (int i{0}; i < k; i++) {
      for (int i1{0}; i1 < loop_ub; i1++) {
        Q[i1 + Q.size(0) * i] = b_A[i1 + nx * i];
      }
    }
  } else {
    binary_expand_op_1(Q, q, s);
  }
  if (r.size(0) == s.size(0)) {
    loop_ub = r.size(0);
    nx = r.size(1);
    b_A.set_size(r.size(0), r.size(1));
    for (int i{0}; i < nx; i++) {
      scalarLB = (r.size(0) / 2) << 1;
      vectorUB = scalarLB - 2;
      for (int i1{0}; i1 <= vectorUB; i1 += 2) {
        __m128d r1;
        b_r = _mm_loadu_pd(&r[i1 + r.size(0) * i]);
        r1 = _mm_loadu_pd(&s[i1]);
        _mm_storeu_pd(&b_A[i1 + b_A.size(0) * i], _mm_mul_pd(b_r, r1));
      }
      for (int i1{scalarLB}; i1 < loop_ub; i1++) {
        b_A[i1 + b_A.size(0) * i] = r[i1 + r.size(0) * i] * s[i1];
      }
    }
    loop_ub = b_A.size(0);
    nx = b_A.size(1);
    r.set_size(b_A.size(0), b_A.size(1));
    for (int i{0}; i < nx; i++) {
      for (int i1{0}; i1 < loop_ub; i1++) {
        r[i1 + r.size(0) * i] = b_A[i1 + b_A.size(0) * i];
      }
    }
    nx = R.size(0);
    for (int i{0}; i < b_loop_ub; i++) {
      for (int i1{0}; i1 < nx; i1++) {
        R[i1 + R.size(0) * i] = r[i1 + nx * i];
      }
    }
  } else {
    binary_expand_op(R, r, s);
  }
}

// End of code generation (qr_unique.cpp)
