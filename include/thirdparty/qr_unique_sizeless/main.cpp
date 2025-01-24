//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// main.cpp
//
// Code generation for function 'main'
//

/*************************************************************************/
/* This automatically generated example C++ main file shows how to call  */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

// Include files
#include "main.h"
#include "coder_array.h"
#include "qr_unique.h"
#include "qr_unique_terminate.h"
#include "rt_nonfinite.h"

#include <iostream>

// Function Declarations
static coder::array<double, 2U>
argInit_UnboundedxUnbounded_real_T(const Eigen::MatrixXd &m);

static double argInit_real_T(const Eigen::MatrixXd &m, int rowId, int colId);

// int main(int argc, char **argv)
// {
//   Eigen::MatrixXd Aeig(4, 3);
//   Aeig(0, 0) = -0.6976;
//   Aeig(1, 0) = -0.7769;
//   Aeig(2, 0) = 0.5055;
//   Aeig(3, 0) = 0.3134;
//   Aeig(0, 1) = -0.5185;
//   Aeig(1, 1) = 0.8990;
//   Aeig(2, 1) = 0.6155;
//   Aeig(3, 1) = 0.1723;
//   Aeig(0, 2) = 0.3196;
//   Aeig(1, 2) = 0.1957;
//   Aeig(2, 2) = -0.1553;
//   Aeig(3, 2) = 0.9518;

//   Eigen::MatrixXd Qeig, Reig;

//   // The initialize function is being called automatically from your entry-point
//   // function. So, a call to initialize is not included here. Invoke the
//   // entry-point functions.
//   // You can call entry-point functions multiple times.
//   main_qr_unique(Aeig, Qeig, Reig);
//   // Terminate the application.
//   // You do not need to do this more than one time.

//   std::cout << "Reig" << std::endl << Reig << std::endl;
//   std::cout << "Qeig" << std::endl << Qeig << std::endl;


//   qr_unique_terminate();
//   return 0;
// }

void main_qr_unique(const Eigen::MatrixXd &Aeig, Eigen::MatrixXd &Qeig,
                    Eigen::MatrixXd &Reig)
{

  std::cout << "Running main_qr_unique" << std::endl;

  coder::array<double, 2U> A;
  coder::array<double, 2U> Q; // Y
  coder::array<double, 2U> R;
  // Initialize function 'qr_unique' input arguments.
  // Initialize function input argument 'A'.
  A = argInit_UnboundedxUnbounded_real_T(Aeig);
  // Call the entry-point 'qr_unique'.
  qr_unique(A, Q, R);

  Qeig.resize(Q.size(0), Q.size(1));
  Reig.resize(R.size(0), R.size(1));

  // std::cout << "R" << std::endl;
  for (int i = 0; i < R.size(0); ++i) {
    for (int j = 0; j < R.size(1); ++j) {
      // std::cout << R[i + R.size(0) * j] << "\t";
      Reig(i,j) = R[i + R.size(0) * j];
    }
    // std::cout << std::endl;
  }

  // std::cout << "Q" << std::endl;
  for (int i = 0; i < Q.size(0); ++i) {
    for (int j = 0; j < Q.size(1); ++j) {
      // std::cout << Q[i + Q.size(0) * j] << "\t";
      Qeig(i,j) = Q[i + Q.size(0) * j];
    }
    // std::cout << std::endl;
  }

  // std::cout << std::endl;
}

// End of code generation (main.cpp)

// Function Definitions
static coder::array<double, 2U>
argInit_UnboundedxUnbounded_real_T(const Eigen::MatrixXd &m)
{
  coder::array<double, 2U> result;
  // Set the size of the array.
  // Change this size to the value that the application requires.
  result.set_size(m.rows(), m.cols());
  // Loop over the array to initialize each element.
  for (int idx0{0}; idx0 < result.size(0); idx0++) {
    for (int idx1{0}; idx1 < result.size(1); idx1++) {
      // Set the value of the array element.
      // Change this value to the value that the application requires.
      result[idx0 + result.size(0) * idx1] = argInit_real_T(m, idx0, idx1);
    }
  }
  return result;
}

static double argInit_real_T(const Eigen::MatrixXd &m, int rowId, int colId)
{
  return m(rowId, colId);
}