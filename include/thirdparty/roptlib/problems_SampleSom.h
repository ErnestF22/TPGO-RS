// #include "manifolds_Rotations.h"
#include "manifolds_Stiefel.h"
#include "manifolds_Euclidean.h"
#include "manifolds_MultiManifolds.h"

#include <vector>
#include <queue>
#include <numeric>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/SparseQR>
#include <eigen3/Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>

#include <rofl/common/io.h>
#include <rofl/common/macros.h>

// #include <ars3d/RshRotation.h>

// #include <ars3d_test/rsh_utils.h>

#include "solvers_RNewton.h" //for RS linesearch

#include "problems_Problem.h"
#include "others_def.h"

#include "manifolds_Manifold.h"
#include "problems_Problem.h"
#include "solvers_RTRNewton.h"

#include "thirdparty/qr_unique_sizeless/main.h"

#include "../../som_utils.h"

namespace ROPTLIB
{

   class SampleSomProblem : public Problem
   {
   public:
      /**
       * Default (empty) constructor
       */
      SampleSomProblem();

      /**
       * Standard constructor with problem data inputs (most commonly used)
       */
      SampleSomProblem(SomUtils::SomSize somSz, SomUtils::MatD &Tijs, Eigen::MatrixXi &edges);

      /**
       * Standard constructor with problem data inputs (most commonly used)
       */
      SampleSomProblem(const SomUtils::SomSize somSz, const SomUtils::MatD &Tijs, const Eigen::MatrixXi &edges);

      /**
       * Standard (empty) destructor
       */
      virtual ~SampleSomProblem();

      /**
       * Compute cost function, given ROPTLIB input x
       * OBS. ROPT to Eig conversion performed here internally
       */
      virtual realdp f(const Variable &x) const;

      /**
       * @brief Compute and return cost with Eigen vectorized input on SE(d)^N
       * (i.e., (d * d * n + d * n) x 1) Eigen input
       */
      double costEigenVecSEdN(const SomUtils::MatD &xEigen) const;

      /**
       * @brief Compute and return cost (as double) with Eigen inputs
       * i.e., p x d x n @param Reigen and p x n @param Teigen
       */
      double costEigen(const SomUtils::VecMatD &Reigen, const SomUtils::MatD &Teigen) const;

      /**
       * @brief Compute and return cost (as double) with fully vectorized
       * (i.e., (p * d * n + p * n) x 1) @param xEigen input
       */
      double costEigenVec(const SomUtils::MatD &xEigen) const;

      /**
       * Computation of Euclidean Gradient of cost function
       */
      // virtual Vector &EucGrad(const Variable &x, Vector *result) const;

      /**
       * Riemannian Gradient (backup: not used, as ROPTLIB wants Grad() directly)
       */
      virtual Vector &RieGrad(const Variable &x, Vector *result) const;

      /**
       * Computation of Riemannian Gradient of cost function
       */
      // virtual Vector &Grad(const Variable &x, Vector *result) const;

      /**
       * Function that computes Euclidean gradient of Translation estimation cost (with Eigen inputs/outputs)
       */
      void egradR(const SomUtils::MatD &P, SomUtils::VecMatD &egR) const;

      /**
       * Function that computes Riemannian gradient of Rotation estimation cost (with Eigen inputs/outputs)
       */
      void rgradR(const SomUtils::VecMatD &R, const SomUtils::MatD &P, SomUtils::VecMatD &rgR) const;

      /**
       * Function that computes Euclidean gradient of Translation estimation cost (with Eigen inputs/outputs)
       */
      void egradT(const SomUtils::MatD &T, const SomUtils::MatD &Lr, const SomUtils::MatD &Pr, SomUtils::MatD &egT) const;

      /**
       * Function that computes Riemannian gradient of Translation estimation cost (with Eigen inputs/outputs)
       */
      void rgradT(const SomUtils::MatD &T, const SomUtils::MatD &Lr, const SomUtils::MatD &Pr, SomUtils::MatD &egT) const;

      /**
       * Euclidean Hessian action i.e., *result = H(x)[etax]
       */
      // virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

      /**
       * Hessian (Genproc) with Eigen I/O;
       * called internally by RieHessianEta()
       */
      void hessGenprocEigen(const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR, const SomUtils::MatD &xT, const SomUtils::MatD &uT, SomUtils::VecMatD &rhR, SomUtils::MatD &rhT) const;

      /**
       * Hessian (Genproc) with Eigen I/O of f-(mu*u)
       */
      void hessGenprocEigenShifted(const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR, const SomUtils::MatD &xT, const SomUtils::MatD &uT, double mu, SomUtils::VecMatD &rhR, SomUtils::MatD &rhT) const;

      /**
       * Hessian action i.e., *result = H(x)[etax]
       */
      virtual Vector &RieHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

      /**
       * Convert ROPTLIB Vector into Eigen equivalent
       * TOCHECK: for now seems to work only with vectors (not matrices, tensors)
       */
      void RoptToEig(Vector x, SomUtils::MatD &xEigen) const;

      /**
       * Get i-th Rotation from ROPTLIB variable x
       */
      void getRi(const Variable &x, SomUtils::MatD &rOut, int i) const;

      /**
       * Get i-th Rotation from Eigen-converted variable xEig
       */
      void getRi(const SomUtils::MatD &xEig, SomUtils::MatD &rOut, int i) const;

      /**
       * Get i-th Rotation from Eigen-converted variable xEig on SE(d)^N
       */
      void getRiSEdN(const SomUtils::MatD &xEig, SomUtils::MatD &rOut, int i) const;

      /**
       * Get all roatations in vector of pxd matrices tOut (with n elements)
       */
      void getRotations(const SomUtils::MatD &xEig, SomUtils::VecMatD &rOut) const;

      /**
       * Get i-th Translation from ROPTLIB variable x
       */
      void getTi(const Variable &x, SomUtils::MatD &rOut, int i) const;

      /**
       * Get i-th Translation from Eigen-converted variable xEig
       */
      void getTi(const SomUtils::MatD &xEig, SomUtils::MatD &tOut, int i) const;

      /**
       * Get i-th Translation from Eigen-converted variable xEig on SE(d)^N
       */
      void getTiSEdN(const SomUtils::MatD &xEig, SomUtils::MatD &tOut, int i) const;

      /**
       * Get all translations in p * n matrix tOut
       */
      void getTranslations(const SomUtils::MatD &xEig, SomUtils::MatD &tOut) const;

      /**
       * Computes matrices used in rotation estimation cost
       */
      void makePfrct(const SomUtils::MatD &T, SomUtils::MatD &P, double &frct) const;

      /**
       * Computes matrices used in translation estimation cost
       */
      void makeLrPrBr(const SomUtils::VecMatD &R, SomUtils::MatD &Lr, SomUtils::MatD &Pr, SomUtils::MatD &Br) const;

      /**
       * Compute one of the Genproc Hessian subparts
       */
      void computeHrr(const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR, const SomUtils::MatD &P, SomUtils::VecMatD &hRR) const;

      /**
       * Compute one of the Genproc Hessian subparts
       */
      void computeHtt(const SomUtils::MatD &uT, const SomUtils::MatD &LR, SomUtils::MatD &hTT) const;

      /**
       * Compute one of the Genproc Hessian subparts
       */
      void computeHrt(const SomUtils::VecMatD &xR, const SomUtils::MatD uT, SomUtils::VecMatD &hrt) const;

      /**
       * Compute one of the Genproc Hessian subparts
       */
      void computeHtr(const SomUtils::VecMatD &uR, SomUtils::MatD &htr) const;

      /**
       * Return p * d (size of a sigle Stiefel-rotation)
       */
      int getRotSz() const;

      /**
       * Return p (size of a sigle Stiefel-translation)
       */
      int getTranslSz() const;

      /**
       * @brief Return vectorized (following col-major order, like in MATLAB) version of @param R
       * in output reference @param RvecOut
       */
      void vectorizeR(const SomUtils::VecMatD &R, SomUtils::MatD &RvecOut) const;

      /**
       * @brief Return vectorized (following col-major order, like in MATLAB) version of @param R
       * Followed by vectorized version of @param T
       * in output reference @param XvecOut
       */
      void vectorizeRT(const SomUtils::VecMatD &R, const SomUtils::MatD &T, SomUtils::MatD &XvecOut) const;

      /**
       * @brief Convert ROPTLIB Vector input @param x on Stiefel manifold into output reference @param xEigen
       */
      void RoptToEigStiefel(Vector x, SomUtils::MatD &xEigen) const;

      /**
       * @brief Make adjacency matrix from class members edges_, sz_.n_
       * and return it in reference @param adjMat
       */
      void makeAdjMatFromEdges(Eigen::MatrixXi &adjMat) const;

      /**
       * @brief Set the Rgt_ (i.e., rotation ground truth) object equal to input @param R
       */
      void setGtR(const SomUtils::VecMatD &R);

      /**
       * @brief Set the Tgt_ (i.e., translation ground truth) object equal to input @param T
       */
      void setGtT(const SomUtils::MatD &T);

      /**
       * @brief Set for the Rgt_ (i.e., rotation ground truth) and Tgt_ (i.e., translation ground truth)
       * objects equal to inputs @param R, @param T
       */
      void setGt(const SomUtils::VecMatD &R, const SomUtils::MatD &T);

      /**
       * @brief Return the Riemannian Staircase Recovery Success (called rsRecoverySuccess_) object
       */
      bool getRsRecoverySuccess() const;

      /**
       * @brief Set the costCurr_ object to @param cc
       */
      void setCostCurr(double cc);

      // private: //TODO: separate public from private members

      ////////////////////////////////////CLASS MEMBER VARIABLES////////////////////////////////////////////////

      /**
       * Struct that contains problem size info
       */
      SomUtils::SomSize sz_;

      /**
       * (Estimated) Relative translations between nodes
       * Size: d x e
       */
      SomUtils::MatD Tijs_;

      /**
       * Edges
       * Size: e x 2
       */
      Eigen::MatrixXi edges_;

      /**
       * Basically just edges_.rows()
       * Saved separately for ease
       */
      int numEdges_; // num edges

      /**
       * Full size of problem: pxdxn (rotations) + pxn (translations)
       */
      int fullSz_;

      /**
       * Ground truth information for R
       */
      SomUtils::VecMatD Rgt_;

      /**
       * Ground truth information for T
       */
      SomUtils::MatD Tgt_;

      /*
       * Boolean stating whether recovery from RS output back to SE(d)^N has been successful or not
       */
      bool rsRecoverySuccess_;

      /**
       * "Global" reference node id (generally 0)
       */
      int src_;

      /**
       * @brief Current cost value (useful at some points e.g. linesearches)
       */
      double costCurr_;

      /**
       * Output of RSOM for R
       */
      SomUtils::VecMatD Rout_;

      /**
       * Output of RSOM for T
       */
      SomUtils::MatD Tout_;

      ////////////////////////////////////////RS////////////////////////////////////////

      /**
       * @brief Check and return boolean stating whether @param lambda and u are an eigencouple for
       * Hessian H(x)[u]
       * where
       * x is composed of rotation and translation components @param xR, @param xT
       * u is composed of rotation and translation components @param uR, @param uT
       * are an eigencouple for hessian at x (up to thr)
       * @return true if it actually is an eigencouple, @return false otherwise
       */
      bool eigencheckHessianGenproc(const double &lambda,
                                    const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR,
                                    const SomUtils::MatD &xT, const SomUtils::MatD &uT, double thr = 1e-3) const;

      /**
       * @brief Check and return boolean stating whether @param lambda and u are an eigencouple for
       * shifted Hessian H(x)[u] - mu*eye()
       * where
       * x is composed of rotation and translation components @param xR, @param xT
       * u is composed of rotation and translation components @param uR, @param uT
       * are an eigencouple for hessian at x (up to thr)
       * @return true if it actually is an eigencouple, @return false otherwise
       */
      bool eigencheckHessianGenprocShifted(const double &lambda,
                                           const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR,
                                           const SomUtils::MatD &xT, const SomUtils::MatD &uT, double mu, double thr = 1e-3) const;
      /**
       * @brief Check whether @param mTg is in the tangent space of @param m
       */
      bool checkIsStiefelTg(const SomUtils::MatD &m, const SomUtils::MatD &mTg) const;

      /**
       * @brief Check whether @param mTg is in the tangent space of @param m
       * 3D Stiefel version
       */
      bool checkIsStiefelTg(const SomUtils::VecMatD &m, const SomUtils::VecMatD &mTg) const;

      /**
       * @brief Check whether @param mTg is in the tangent space of @param m and has norm 1
       */
      bool checkIsStiefelTgNorm(const SomUtils::MatD &m, const SomUtils::MatD &mTg) const;

      /**
       * @brief Generate a random normalized vector tangent @param mIn on 2D Stiefel manifold
       * and return it in reference @param mOut
       */
      void stiefelRandTgNormVector(const SomUtils::MatD &mIn, SomUtils::MatD &mOut) const;

      /**
       * @brief Generate a random normalized vector tangent @param mIn on 3D Stiefel manifold
       * and return it in reference @param mOut
       */
      void stiefelRandTgNormVector(const SomUtils::VecMatD &mIn, SomUtils::VecMatD &mOut) const;

      /**
       * @brief Power iteration method of the genproc hessian
       * Hessian is considered as H(x)[u] where
       * x is composed of rotation and translation components @param xR, @param xT
       * u is composed of rotation and translation components @param uR, @param uT
       * @param lambdaMax @param uOutR @param uOutT are the reference outputs for
       * @param thresh is used for early stopping conditions
       * Number of max iterations of PIM is set empirically to 2500
       */
      void pimFunctionGenproc(const SomUtils::VecMatD &xR, const SomUtils::MatD &xT, const SomUtils::VecMatD &uR, const SomUtils::MatD &uT,
                              double &lambdaMax, SomUtils::VecMatD &uOutR, SomUtils::MatD &uOutT, double thresh = 1e-5) const;

      /**
       * @brief Power iteration method of the genproc hessian shifted by @param mu
       * Hessian is considered as: H(x)[u] - mu*eye()
       * where
       * x is composed of rotation and translation components @param xR, @param xT
       * u is composed of rotation and translation components @param uR, @param uT
       * @param lambdaMax @param uOutR @param uOutT are the reference outputs for
       * @param thresh is used for early stopping conditions
       * Number of max iterations of PIM is set empirically to 2500
       */
      void pimFunctionGenprocShifted(const SomUtils::VecMatD &xR, const SomUtils::MatD &xT, const SomUtils::VecMatD &uR, const SomUtils::MatD &uT, double mu,
                                     double &lambdaMax, SomUtils::VecMatD &uOutR, SomUtils::MatD &uOutT, double thresh = 1e-5) const;

      /**
       * For debug purposes, reduced version of rsomPimHessianGenproc, up to linesearch (excluded)
       * No outputs, just log prints
       *
       * @brief Power iteration method for R, T (generalized Procrustes) version of the problem
       * Params are:
       * @param thresh for PIM computation thresholds (eigencouple check, stopping conditions of PIM)
       * @param R, @param T are the starting points
       */
      void rsomPimHessianGenprocSmall(double thresh, const SomUtils::VecMatD &R, const SomUtils::MatD &T) const;

      /**
       * @brief Power iteration method for R, T (generalized Procrustes) version of the problem
       * Params are:
       * @param thresh for PIM computation thresholds (eigencouple check, stopping conditions of PIM)
       * @param R, @param T are the starting points
       * @param Y0 reference output: new starting point for RS next step
       * @param lambdaPimOut, vPimRout, vPimTout reference outputs: associated eigencouple
       * @param armijo (= false by default) indicates whether the linesearch method used is based
       * on ROPTLIB's Armijo-Goldstein (armijo = true) implementation of simply linesearchDummy (armijo = false)
       */
      void rsomPimHessianGenprocEigen(double thresh,
                                      const SomUtils::VecMatD &R, const SomUtils::MatD &T,
                                      Vector &Y0, double &lambdaPimOut, SomUtils::VecMatD &vPimRout, SomUtils::MatD &vPimTout,
                                      bool armijo = false) const;

      /**
       * @brief Find minimum eigenvalue of Hessian H(x)[u] using a basis of the tangent space
       * @param R, @param T are the starting points
       * @param Y0 reference output: new starting point for RS next step
       * @param lambdaPimOut, vPimRout, vPimTout reference outputs: associated eigencouple
       */
      void rsomEscapeHessianGenprocEigen(const SomUtils::VecMatD &R, const SomUtils::MatD &T,
                                         Vector &Y0, double &lambdaPimOut, SomUtils::VecMatD &vPimRout, SomUtils::MatD &vPimTout) const;

      /**
       * @brief Power iteration method for R, T (generalized Procrustes) version of the problem
       * Params are:
       * @param thresh for PIM computation thresholds (eigencouple check, stopping conditions of PIM)
       * @param R, @param T are the starting points
       * @param Y0 reference output
       * @param armijo (= false by default) indicates whether the linesearch method used is based
       * on ROPTLIB's Armijo-Goldstein (armijo = true) implementation of simply linesearchDummy (armijo = false)
       */
      void rsomPimHessianGenproc(double thresh, const SomUtils::VecMatD &R, const SomUtils::MatD &T, Vector &Y0, bool armijo = false) const;

      /**
       * @brief Linesearch for cost decrease based on Armijo conditions (see ROPTLIB documentation) from @param xIn starting point
       * with size @param somSzLocal returning output in reference @param Y0
       */
      void linesearchArmijoROPTLIB(const Vector &xIn, const SomUtils::SomSize &somSzLocal, Vector &Y0) const;

      /**
       * @brief Thin QR factorization ensuring diagonal of R is real, positive if possible.
          % If A is a matrix, then Q, R are matrices such that A = QR where Q'*Q = I
          % and R is upper triangular. If A is real, then so are Q and R.
          % This is a thin QR factorization in the sense that if A has more rows than
          % columns, than Q has the same size as A.
          %
          % If A has full column rank, then R has positive reals on its diagonal.
          % Otherwise, it may have zeros on its diagonal.
          %
          % This is equivalent to a call to Matlab's qr(A, 0), up to possible
          % sign/phase changes of the columns of Q and of the rows of R to ensure
          % the stated properties of the diagonal of R. If A has full column rank,
          % this decomposition is unique, hence the name of the function.
          %
          % If A is a 3D array, then Q, R are also 3D arrays and this function is
          % applied to each slice separately.
          % Adapted from Manopt: www.manopt.org.

          2D version
       */
      //   void QRunique(const SomUtils::MatD &A, SomUtils::MatD &Q, SomUtils::MatD &R) const;

      /**
       * @brief Thin QR factorization ensuring diagonal of R is real, positive if possible.
          % If A is a matrix, then Q, R are matrices such that A = QR where Q'*Q = I
          % and R is upper triangular. If A is real, then so are Q and R.
          % This is a thin QR factorization in the sense that if A has more rows than
          % columns, than Q has the same size as A.
          %
          % If A has full column rank, then R has positive reals on its diagonal.
          % Otherwise, it may have zeros on its diagonal.
          %
          % This is equivalent to a call to Matlab's qr(A, 0), up to possible
          % sign/phase changes of the columns of Q and of the rows of R to ensure
          % the stated properties of the diagonal of R. If A has full column rank,
          % this decomposition is unique, hence the name of the function.
          %
          % If A is a 3D array, then Q, R are also 3D arrays and this function is
          % applied to each slice separately.
          % Adapted from Manopt: www.manopt.org.

          3D version
       */
      void QRunique(const SomUtils::VecMatD &A, SomUtils::VecMatD &Q, SomUtils::VecMatD &R) const;

      /**
       * @brief Return whether input @param m is on 3D Stiefel manifold
       */
      bool checkIsOn3dStiefel(const SomUtils::VecMatD &m) const;

      /**
       * @brief Return whether input @param m is on 2D Stiefel manifold
       */
      bool checkIsOnStiefel(const SomUtils::MatD &m) const;

      /**
       * @brief Stiefel Retraction for 3D manifold from @param xIn towards direction @param e multiplied by @param t
       * returning new point in reference @param rxe
       * Retraction version based on QR "unique" decomposition
       */
      void stiefelRetractionQR(const SomUtils::VecMatD &x, const SomUtils::VecMatD &e, SomUtils::VecMatD &rxe, double t = 1.0) const;

      /**
       * @brief Euclidean Retraction for 2D manifold from @param x towards direction @param d with @param t steps (default t = 1.0)
       * returning new point in reference @param y
       * Retraction version based on sqrtm
       */
      void euclRetraction(const SomUtils::MatD &x, const SomUtils::MatD &d, SomUtils::MatD &y, double t = 1.0) const;

      /**
       * @brief Stiefel Retraction for 2D manifold from @param xIn towards direction @param e
       * returning new point in reference @param rxe
       * Retraction version based on sqrtm
       */
      void stiefelRetractionPolar(const SomUtils::MatD &xIn, const SomUtils::MatD &e, SomUtils::MatD &rxe, double t = 1.0) const;

      /**
       * @brief Stiefel Retraction for 3D manifold from @param xIn towards direction @param e
       * returning new point in reference @param rxe
       */
      void stiefelRetractionPolar(const SomUtils::VecMatD &xIn, const SomUtils::VecMatD &e, SomUtils::VecMatD &rxe, double t = 1.0) const;

      /**
       * @brief Call stiefelRetractionPolar, euclRetraction on szNext elements @param xRin, @param xTin, @param vRin, @param vTin
       * Return new point "Y0" in reference @param Y0R and @param Y0T
       */
      void linesearchDummy(const double costInit,
                           const SomUtils::VecMatD &xRin, const SomUtils::MatD &xTin,
                           const SomUtils::VecMatD &vRin, const SomUtils::MatD &vTin,
                           SomUtils::VecMatD &Y0R, SomUtils::MatD &Y0T,
                           bool qr = true) const;

      ////////////////////////////////////////RECOVERY////////////////////////////////////////

      /**
       * @brief Compute (and return them in reference vector @param nodeDegrees)
       * node degrees from edges_ class member (adjacency matrix is built internally)
       */
      void computeNodeDegrees(Eigen::ArrayXi &nodeDegrees) const;

      /**
       * @brief Return in reference output @param Tedges
       * an array with differences T(:,i) - T(:,j) for all (i,j) edge pairs in edges_
       * This version of the function also returns @param T1offset = T.col(src_);
       */
      void makeTedges(const SomUtils::MatD &T, SomUtils::MatD &Tedges, SomUtils::MatD &T1offset) const;

      /**
       * @brief Return in reference output @param Tedges
       * an array with differences T(:,i) - T(:,j) for all (i,j) edge pairs in edges_
       */
      void makeTedges(const SomUtils::MatD &T, SomUtils::MatD &Tedges) const;

      /**
       * @brief Return in reference @param Qtransp an aligned (i.e., with last row to 0) version of input @param x
       * through svd (if possible)
       */
      void POCRotateToMinimizeLastEntries(const SomUtils::MatD &x, SomUtils::MatD &Qtransp) const;

      /**
       * @brief Make Tij1j2, \tilde{Tij1j2} and return them to reference params @param Tij1j2, @param Tij1j2tilde
       * for node @param nodeId
       * using input params @param nodeDegrees, @param Tedges accordingly
       */
      void makeTij1j2sEdges(int nodeId, const Eigen::ArrayXi &nodeDegrees, const SomUtils::MatD &Tedges,
                            SomUtils::MatD &Tij1j2, SomUtils::MatD &Tij1j2tilde) const;

      /**
       * @brief Compute a basis for the orthogonal complement to the columns of @param v, i.e.,
       * vOrth'*v=0 and return result to reference @param vOrth
       */
      void orthComplement(const SomUtils::MatD &v, SomUtils::MatD &vOrth) const;
      /*
       * @brief Complete the columns of @param Q (assumed orthonormal) to an orthonormal basis
       * (output is reference @param Qout)
       */
      void orthCompleteBasis(const SomUtils::MatD &Q, SomUtils::MatD &Qout) const;

      /**
       * @brief Flip input vector @param Ain from left-to-right and return result to reference @param Aout
       */
      void fliplr(const SomUtils::MatD &Ain, SomUtils::MatD &Aout) const;

      /**
       * @brief Flip input vector @param Ain upside-down and return result to reference @param Aout
       */
      void flipud(const SomUtils::MatD &Ain, SomUtils::MatD &Aout) const;

      /**
       * @brief Util function used internally by recoverRiTilde
       */
      void align2d(const SomUtils::MatD &v, SomUtils::MatD &Qx) const;

      /**
       * @brief Util function used internally by recoverRiTilde
       */
      void procrustesRb(const SomUtils::MatD &c, const SomUtils::MatD &q, SomUtils::MatD &RbEst) const;

      /**
       * @brief Recover a SO(d) solution from RS output @param RiTilde2
       * using @param Tijtilde (corresponding elements in Tijs_)
       *
       * This algorithm returns two possible solutions
       * @param RiTildeEst1, @param RiTildeEst2
       * which can be discerned only later when solving gauge symmetries
       */
      void recoverRiTilde(const SomUtils::MatD &RiTilde2, const SomUtils::MatD &Tijtilde,
                          SomUtils::MatD &RiTildeEst1, SomUtils::MatD &RiTildeEst2) const;

      /**
       * @brief Dummy (almost brute-force) version of the backtracking algorithm that, after calling dijkstraSP()
       * gives back in output referece @param listEdges
       * a list of nodes (identified through ids) linking source node @param src to the @param n nodes of the graph
       * via a shortest path
       *
       * @param prev is output of dijkstraSP()
       */
      void dijkstraBT(int src, int n, const std::vector<int> &prev, std::vector<std::vector<int>> &list) const;

      /**
       * @brief Dummy (almost brute-force) version of the backtracking algorithm that, after calling dijkstraSP()
       * gives back in output referece @param listEdges
       * a list of edges (identified through ids) linking source node @param src to the @param n nodes of the graph
       * via a shortest path
       *
       * @param prev, @param edges are outputs of dijkstraSP()
       */
      void dijkstraBTedges(int src, int n, const std::vector<int> &prev, const Eigen::MatrixXi &edges, std::vector<std::vector<int>> &listEdges) const;

      /**
       * @brief Implementation of the Dijkstra classic shortest path algorithm on a graph defined by adjacency matrix @param adjmat
       * with @param n nodes, setting source node to the one with id = @param src as reference node
       * Usual outputs (distance from source, previous node in path per each node)
       * are returned in reference @param dist and reference @param prev
       */
      void dijkstraSP(int n, int src, const Eigen::MatrixXi &adjmat, std::vector<double> &dist, std::vector<int> &prev) const;

      /**
       * @brief From edge differences vector @param Tdiffs
       * build output reference vector @param T composed of @param n nodes/elements
       * with global reference node having id = @param src
       */
      void edgeDiffs2T(int src, const SomUtils::MatD &Tdiffs, int n, SomUtils::MatD &T) const;

      /**
       * @brief Recover RS output composed of @param RmanoptOut, @param TmanoptOut
       * after RS exit at @param staircaseStepIdx
       * Return output is reference @param Rrecovered and @param Trecovered
       * @return true if recovery is successful (determinant check); @return false otherwise
       */
      bool recoverySEdN(int staircaseStepIdx,
                        const SomUtils::VecMatD &RmanoptOut, const SomUtils::MatD &TmanoptOut,
                        SomUtils::VecMatD &Rrecovered, SomUtils::MatD &Trecovered);

      /**
       * @brief Solve global gauge for recoverySEdN() outputs @param Rsedn and @param Tsedn,
       * setting @param src node rotation and translation to a certain known quantity (ground truth info)
       * and rototranslate all the other rotations and translations accordingly in order to maintain
       * reciprocal distances between nodes
       * @return true if "globalization" is successful (determinant check); @return false otherwise
       */
      bool globalize(int src, const SomUtils::VecMatD &Rsedn, const SomUtils::MatD &Tsedn,
                     SomUtils::VecMatD &Rout, SomUtils::MatD &Tout);

      /**
       * @brief Make Hmat matrix s.t. Hmat * u = H(x)[u] for all u in tangent space
       * of x = [R, T] (with R in SO(d)^N, T in R^(d x N))
       * @param XvecNext is the vectorized version of x
       * @param szNext is the size of the problem
       * @param Tijs is the relative translations between nodes
       * @param edges is the adjacency matrix
       * reference @param Hmat is the output matrix
       */
      void makeHmat(const SomUtils::MatD &XvecNext, const SomUtils::SomSize &szNext, SomUtils::MatD &Hmat) const;
   };

   /**
    * @brief Solve an instance of the RSOM problem using ICP-based method
    */
   double runRsomICP(ROPTLIB::SampleSomProblem &Prob, const ROPTLIB::Vector &startX,
                     int src, SomUtils::VecMatD &Rout, SomUtils::MatD &Tout,
                     int numMaxIter = 20, double stopThr = 1e-3);

   /**
    * @brief Solve an instance of the RSOM problem using RS
    */
   double runRsomRS(ROPTLIB::SampleSomProblem &Prob, const ROPTLIB::Vector &startX,
                    int src,
                    SomUtils::VecMatD &Rout, SomUtils::MatD &Tout,
                    int &staircaseStepIdx);

} // end of namespace ROPTLIB

// "CODE FOR PLOTTING CONCAVITY"
// alphas = linspace(-0.01,0.01,501); %-0.2:0.01:0.2;
// plot_vals = zeros(size(alphas));
// plot_vals_taylor = zeros(size(alphas));
// for ii = 1:length(alphas)
//     x_retr_ii = step2.M.retr(Xnext, v_pim_after_shift, alphas(ii));
// %     disp("Is x_retr_ii on Stiefel? (Taylor)")
// %     disp(check_is_on_stiefel(x_retr_ii));
// %     disp([matStack(x), matStack(x_retr_ii)])
//     plot_vals(ii) = step2.cost(x_retr_ii);
//     %Note: gradient is zero
//     pvt_R = alphas(ii)^2/2* ...
//         sum(stiefel_metric( ...
//         Rnext,v_pim_after_shift.R, ...
//             hess_genproc_R(Xnext, v_pim_after_shift, problem_struct_next),'canonical'));
//     pvt_T = alphas(ii)^2/2* ...
//         sum(stiefel_metric( ...
//         Tnext,v_pim_after_shift.T, ...
//             hess_genproc_T(Xnext, v_pim_after_shift, problem_struct_next),'euclidean'));
//     plot_vals_taylor(ii) = step2.cost(Xnext) + pvt_R + pvt_T;
// end

// plot(alphas, plot_vals,'b')
// hold on
// plot(alphas,plot_vals_taylor,'k.');
// hold off