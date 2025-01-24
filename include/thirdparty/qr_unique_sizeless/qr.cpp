//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// qr.cpp
//
// Code generation for function 'qr'
//

// Include files
#include "qr.h"
#include "rt_nonfinite.h"
#include "xnrm2.h"
#include "coder_array.h"
#include <cmath>
#include <emmintrin.h>

// Function Definitions
namespace coder {
void qr(const array<double, 2U> &A, array<double, 2U> &Q, array<double, 2U> &R)
{
  __m128d r;
  array<double, 2U> b_A;
  array<double, 1U> tau;
  array<double, 1U> work;
  double beta1;
  double temp;
  int b_i;
  int exitg1;
  int i;
  int i1;
  int ia;
  int ii;
  int itau;
  int jA;
  int knt;
  int lastv;
  int m;
  int minsz;
  int mmi_tmp;
  int vectorUB_tmp;
  boolean_T exitg2;
  i = A.size(0);
  minsz = A.size(1);
  b_A.set_size(A.size(0), A.size(1));
  b_i = A.size(0) * A.size(1);
  for (i1 = 0; i1 < b_i; i1++) {
    b_A[i1] = A[i1];
  }
  m = A.size(0) - 1;
  b_i = A.size(0);
  itau = A.size(1);
  if (b_i <= itau) {
    itau = b_i;
  }
  tau.set_size(itau);
  for (i1 = 0; i1 < itau; i1++) {
    tau[i1] = 0.0;
  }
  if ((A.size(0) != 0) && (A.size(1) != 0) && (itau >= 1)) {
    work.set_size(A.size(1));
    for (i1 = 0; i1 < minsz; i1++) {
      work[i1] = 0.0;
    }
    for (int c_i{0}; c_i < itau; c_i++) {
      double atmp;
      ii = c_i * i + c_i;
      mmi_tmp = m - c_i;
      if (c_i < m) {
        atmp = b_A[ii];
        lastv = ii + 2;
        tau[c_i] = 0.0;
        if (mmi_tmp + 1 > 0) {
          beta1 = internal::blas::xnrm2(mmi_tmp, b_A, ii + 2);
          if (beta1 != 0.0) {
            temp = std::abs(b_A[ii]);
            beta1 = std::abs(beta1);
            if (temp < beta1) {
              temp /= beta1;
              beta1 *= std::sqrt(temp * temp + 1.0);
            } else if (temp > beta1) {
              beta1 /= temp;
              beta1 = temp * std::sqrt(beta1 * beta1 + 1.0);
            } else if (std::isnan(beta1)) {
              beta1 = rtNaN;
            } else {
              beta1 = temp * 1.4142135623730951;
            }
            if (b_A[ii] >= 0.0) {
              beta1 = -beta1;
            }
            if (std::abs(beta1) < 1.0020841800044864E-292) {
              knt = 0;
              i1 = (ii + mmi_tmp) + 1;
              do {
                knt++;
                b_i = (((((i1 - ii) - 1) / 2) << 1) + ii) + 2;
                vectorUB_tmp = b_i - 2;
                for (jA = lastv; jA <= vectorUB_tmp; jA += 2) {
                  r = _mm_loadu_pd(&b_A[jA - 1]);
                  _mm_storeu_pd(
                      &b_A[jA - 1],
                      _mm_mul_pd(_mm_set1_pd(9.9792015476736E+291), r));
                }
                for (jA = b_i; jA <= i1; jA++) {
                  b_A[jA - 1] = 9.9792015476736E+291 * b_A[jA - 1];
                }
                beta1 *= 9.9792015476736E+291;
                atmp *= 9.9792015476736E+291;
              } while ((std::abs(beta1) < 1.0020841800044864E-292) &&
                       (knt < 20));
              temp = std::abs(atmp);
              beta1 = std::abs(internal::blas::xnrm2(mmi_tmp, b_A, ii + 2));
              if (temp < beta1) {
                temp /= beta1;
                beta1 *= std::sqrt(temp * temp + 1.0);
              } else if (temp > beta1) {
                beta1 /= temp;
                beta1 = temp * std::sqrt(beta1 * beta1 + 1.0);
              } else if (std::isnan(beta1)) {
                beta1 = rtNaN;
              } else {
                beta1 = temp * 1.4142135623730951;
              }
              if (atmp >= 0.0) {
                beta1 = -beta1;
              }
              tau[c_i] = (beta1 - atmp) / beta1;
              temp = 1.0 / (atmp - beta1);
              for (jA = lastv; jA <= vectorUB_tmp; jA += 2) {
                r = _mm_loadu_pd(&b_A[jA - 1]);
                _mm_storeu_pd(&b_A[jA - 1], _mm_mul_pd(_mm_set1_pd(temp), r));
              }
              for (jA = b_i; jA <= i1; jA++) {
                b_A[jA - 1] = temp * b_A[jA - 1];
              }
              for (jA = 0; jA < knt; jA++) {
                beta1 *= 1.0020841800044864E-292;
              }
              atmp = beta1;
            } else {
              tau[c_i] = (beta1 - b_A[ii]) / beta1;
              temp = 1.0 / (b_A[ii] - beta1);
              i1 = (ii + mmi_tmp) + 1;
              b_i = (((((i1 - ii) - 1) / 2) << 1) + ii) + 2;
              knt = b_i - 2;
              for (jA = lastv; jA <= knt; jA += 2) {
                r = _mm_loadu_pd(&b_A[jA - 1]);
                _mm_storeu_pd(&b_A[jA - 1], _mm_mul_pd(_mm_set1_pd(temp), r));
              }
              for (jA = b_i; jA <= i1; jA++) {
                b_A[jA - 1] = temp * b_A[jA - 1];
              }
              atmp = beta1;
            }
          }
        }
        b_A[ii] = atmp;
      } else {
        tau[c_i] = 0.0;
      }
      if (c_i + 1 < minsz) {
        atmp = b_A[ii];
        b_A[ii] = 1.0;
        jA = (ii + i) + 1;
        if (tau[c_i] != 0.0) {
          lastv = mmi_tmp + 1;
          b_i = ii + mmi_tmp;
          while ((lastv > 0) && (b_A[b_i] == 0.0)) {
            lastv--;
            b_i--;
          }
          knt = (minsz - c_i) - 1;
          exitg2 = false;
          while ((!exitg2) && (knt > 0)) {
            b_i = jA + (knt - 1) * i;
            ia = b_i;
            do {
              exitg1 = 0;
              if (ia <= (b_i + lastv) - 1) {
                if (b_A[ia - 1] != 0.0) {
                  exitg1 = 1;
                } else {
                  ia++;
                }
              } else {
                knt--;
                exitg1 = 2;
              }
            } while (exitg1 == 0);
            if (exitg1 == 1) {
              exitg2 = true;
            }
          }
        } else {
          lastv = 0;
          knt = 0;
        }
        if (lastv > 0) {
          if (knt != 0) {
            for (vectorUB_tmp = 0; vectorUB_tmp < knt; vectorUB_tmp++) {
              work[vectorUB_tmp] = 0.0;
            }
            vectorUB_tmp = 0;
            i1 = jA + i * (knt - 1);
            for (mmi_tmp = jA; i < 0 ? mmi_tmp >= i1 : mmi_tmp <= i1;
                 mmi_tmp += i) {
              beta1 = 0.0;
              b_i = mmi_tmp + lastv;
              for (ia = mmi_tmp; ia < b_i; ia++) {
                beta1 += b_A[ia - 1] * b_A[(ii + ia) - mmi_tmp];
              }
              work[vectorUB_tmp] = work[vectorUB_tmp] + beta1;
              vectorUB_tmp++;
            }
          }
          beta1 = -tau[c_i];
          if (!(beta1 == 0.0)) {
            for (vectorUB_tmp = 0; vectorUB_tmp < knt; vectorUB_tmp++) {
              if (work[vectorUB_tmp] != 0.0) {
                temp = work[vectorUB_tmp] * beta1;
                i1 = lastv + jA;
                for (b_i = jA; b_i < i1; b_i++) {
                  b_A[b_i - 1] = b_A[b_i - 1] + b_A[(ii + b_i) - jA] * temp;
                }
              }
              jA += i;
            }
          }
        }
        b_A[ii] = atmp;
      }
    }
  }
  m = b_A.size(0);
  knt = b_A.size(1);
  b_i = b_A.size(0);
  minsz = b_A.size(1);
  if (b_i <= minsz) {
    minsz = b_i;
  }
  R.set_size(minsz, b_A.size(1));
  for (vectorUB_tmp = 0; vectorUB_tmp < minsz; vectorUB_tmp++) {
    for (int c_i{0}; c_i <= vectorUB_tmp; c_i++) {
      R[c_i + R.size(0) * vectorUB_tmp] = b_A[c_i + b_A.size(0) * vectorUB_tmp];
    }
    i = vectorUB_tmp + 2;
    for (int c_i{i}; c_i <= minsz; c_i++) {
      R[(c_i + R.size(0) * vectorUB_tmp) - 1] = 0.0;
    }
  }
  i = b_A.size(0) + 1;
  for (vectorUB_tmp = i; vectorUB_tmp <= knt; vectorUB_tmp++) {
    for (int c_i{0}; c_i < minsz; c_i++) {
      R[c_i + R.size(0) * (vectorUB_tmp - 1)] =
          b_A[c_i + b_A.size(0) * (vectorUB_tmp - 1)];
    }
  }
  if (minsz >= 1) {
    for (vectorUB_tmp = minsz; vectorUB_tmp < minsz; vectorUB_tmp++) {
      ia = vectorUB_tmp * m;
      for (int c_i{0}; c_i < m; c_i++) {
        b_A[ia + c_i] = 0.0;
      }
      b_A[ia + vectorUB_tmp] = 1.0;
    }
    itau = minsz - 1;
    work.set_size(b_A.size(1));
    for (i = 0; i < knt; i++) {
      work[i] = 0.0;
    }
    for (int c_i{minsz}; c_i >= 1; c_i--) {
      ii = c_i + (c_i - 1) * m;
      if (c_i < minsz) {
        b_A[ii - 1] = 1.0;
        jA = ii + m;
        if (tau[itau] != 0.0) {
          lastv = (m - c_i) + 1;
          b_i = jA - c_i;
          while ((lastv > 0) && (b_A[b_i - 1] == 0.0)) {
            lastv--;
            b_i--;
          }
          knt = minsz - c_i;
          exitg2 = false;
          while ((!exitg2) && (knt > 0)) {
            b_i = jA + (knt - 1) * m;
            ia = b_i;
            do {
              exitg1 = 0;
              if (ia <= (b_i + lastv) - 1) {
                if (b_A[ia - 1] != 0.0) {
                  exitg1 = 1;
                } else {
                  ia++;
                }
              } else {
                knt--;
                exitg1 = 2;
              }
            } while (exitg1 == 0);
            if (exitg1 == 1) {
              exitg2 = true;
            }
          }
        } else {
          lastv = 0;
          knt = 0;
        }
        if (lastv > 0) {
          if (knt != 0) {
            for (vectorUB_tmp = 0; vectorUB_tmp < knt; vectorUB_tmp++) {
              work[vectorUB_tmp] = 0.0;
            }
            vectorUB_tmp = 0;
            i = jA + m * (knt - 1);
            for (mmi_tmp = jA; m < 0 ? mmi_tmp >= i : mmi_tmp <= i;
                 mmi_tmp += m) {
              beta1 = 0.0;
              i1 = mmi_tmp + lastv;
              for (ia = mmi_tmp; ia < i1; ia++) {
                beta1 += b_A[ia - 1] * b_A[((ii + ia) - mmi_tmp) - 1];
              }
              work[vectorUB_tmp] = work[vectorUB_tmp] + beta1;
              vectorUB_tmp++;
            }
          }
          beta1 = -tau[itau];
          if (!(beta1 == 0.0)) {
            for (vectorUB_tmp = 0; vectorUB_tmp < knt; vectorUB_tmp++) {
              if (work[vectorUB_tmp] != 0.0) {
                temp = work[vectorUB_tmp] * beta1;
                i = lastv + jA;
                for (b_i = jA; b_i < i; b_i++) {
                  b_A[b_i - 1] =
                      b_A[b_i - 1] + b_A[((ii + b_i) - jA) - 1] * temp;
                }
              }
              jA += m;
            }
          }
        }
      }
      if (c_i < m) {
        lastv = ii + 1;
        i = (ii + m) - c_i;
        b_i = ((((i - ii) / 2) << 1) + ii) + 1;
        knt = b_i - 2;
        for (jA = lastv; jA <= knt; jA += 2) {
          r = _mm_loadu_pd(&b_A[jA - 1]);
          _mm_storeu_pd(&b_A[jA - 1], _mm_mul_pd(_mm_set1_pd(-tau[itau]), r));
        }
        for (jA = b_i; jA <= i; jA++) {
          b_A[jA - 1] = -tau[itau] * b_A[jA - 1];
        }
      }
      b_A[ii - 1] = 1.0 - tau[itau];
      for (vectorUB_tmp = 0; vectorUB_tmp <= c_i - 2; vectorUB_tmp++) {
        b_A[(ii - vectorUB_tmp) - 2] = 0.0;
      }
      itau--;
    }
  }
  Q.set_size(b_A.size(0), minsz);
  for (vectorUB_tmp = 0; vectorUB_tmp < minsz; vectorUB_tmp++) {
    for (int c_i{0}; c_i < m; c_i++) {
      Q[c_i + Q.size(0) * vectorUB_tmp] = b_A[c_i + b_A.size(0) * vectorUB_tmp];
    }
  }
}

} // namespace coder

// End of code generation (qr.cpp)
