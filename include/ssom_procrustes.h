#ifndef SSOM_PROCRUSTES_H_
#define SSOM_PROCRUSTES_H_

#include "som_utils.h"
#include "thirdparty/roptlib/cwrapper/lapack/dgesv.h"

class SsomProcrustes
{
public:
    /**
     * Default (empty) constructor
     */
    SsomProcrustes();

    /**
     * Standard constructor with problem data inputs (most commonly used)
     */
    SsomProcrustes(SomUtils::SomSize somSz, SomUtils::MatD &Tijs, Eigen::MatrixXi &edges);

    /**
     * Standard constructor with problem data inputs (most commonly used)
     */
    SsomProcrustes(const SomUtils::SomSize somSz, const SomUtils::MatD &Tijs, const Eigen::MatrixXi &edges);

    /**
     * Default (empty) destructor
     */
    virtual ~SsomProcrustes();

private:
    /**
     * Extract i-th part from Mmat or Nmat (see makeMN)
     */
    void getMiFromMmat(const SomUtils::MatD &Mmat, SomUtils::MatD &Mi, int idx) const;

    /**
     * Estimate rotations
     */
    void stepOne(const SomUtils::MatD &T, const SomUtils::MatD &Lambdas);

    /**
     * Estimate translations
     */
    void stepTwo(const SomUtils::VecMatD &R, const SomUtils::MatD &Lambdas);

    /**
     * Estimate Lambdas
     */
    void stepThree(const SomUtils::VecMatD &R,
                    const SomUtils::MatD &T);

    /**
       * @brief Solve a quadratic equation
       * @param a coefficient a
       * @param b coefficient b
       * @param c coefficient c
       * @param x output root
       */
      double solveQuadraticLambdas(const double a, const double b, const double c);                    

    /**
     * Make matrices M and N (used in stepOne)
     */
    void makeMN(const SomUtils::MatD &T, SomUtils::MatD &Mmat, SomUtils::MatD &Nmat, SomUtils::MatD &TijsScaled) const;

    /**
     * Make scaled Tijs
     */
    void makeTijsScaled(const SomUtils::MatD &Tijs, const SomUtils::MatD &Lambdas, SomUtils::MatD &TijsScaled) const;

public:
    /**
     * Make A and b matrices (used in stepTwo)
     */
    void makeAb(const SomUtils::VecMatD &Rgf, SomUtils::MatD &A, SomUtils::MatD &b, SomUtils::MatD &TijsScaled) const;

    /**
     * Run the entire pipeline
     */
    void run();

    /**
     * Getter for costCurr_ object
     */
    double getCost() const;

    /**
     * @brief Compute and return cost (as double) with Eigen inputs
     * i.e., p x d x n @param Reigen and p x n @param Teigen
     */
    double costEigen(const SomUtils::VecMatD &Reigen, const SomUtils::MatD &Teigen, const SomUtils::MatD &LambdasEigen) const;

    /**
     * Set Tstart_ before running stepOne for the first time
     */
    void setTstart(const SomUtils::MatD &Tstart);

    /**
     * Set LambdasStart_ before running stepThree for the first time
     */
    void setLambdasStart(const SomUtils::MatD &LambdasStart);

    /**
     * Set ground truth information object for R
     */
    void setRgt(const SomUtils::VecMatD &Rgt);

    /**
     * Set ground truth information object for T
     */
    void setTgt(const SomUtils::MatD &Tgt);

    /**
     * Set ground truth information object for Lambdas
     */
    void setLambdasGt(const SomUtils::MatD &LambdasGt);

    /**
     * Get the output of SSOM-Procrustes for T
     */
    SomUtils::MatD getTout() const;

    /**
     * Get the output of SSOM-Procrustes for R
     */
    SomUtils::VecMatD getRout() const;

    /**
     * Get the output of SSOM-Procrustes for Lambdas
     */
    SomUtils::MatD getLambdasOut() const;

private:
    /**
     * Vectorize R, T, and Lambdas into a single vector XvecOut
     * equivalent of Matlab's XvecOut = [R(:); T(:); Lambdas(:)];
     */
    void vectorizeRTLambdas(const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas, SomUtils::MatD &XvecOut) const;

    /**
     * Compute norm of [R1(:); T1(:); Lambdas1(:)] - [R2(:); T2(:); Lambdas2(:)]
     */
    double norm(const SomUtils::VecMatD &R1, const SomUtils::MatD &T1, const SomUtils::MatD &Lambdas1, const SomUtils::VecMatD &R2, const SomUtils::MatD &T2, const SomUtils::MatD &Lambdas2) const;

    /**
     * Compute norm of [R(:); T(:); Lambdas(:)]
     */
    double norm(const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas) const;

    /**
     * Set Tcurr_ object (needed in testing-phase, now private)
     */
    void setTcurr(const SomUtils::MatD &Tcurr);

    /**
     * CLASS MEMBERS
     */

    /**
     * Struct that contains problem size info
     */
    SomUtils::SomSize sz_;

    /**
     * (Estimated) Relative translations between nodes
     * Size: d x e
     */
    SomUtils::MatD tijs_;

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

    /**
     * Ground truth information for Lambdas
     */
    SomUtils::MatD LambdasGt_;

    /**
     * "Global" reference node id (generally 0)
     */
    int src_;

    /**
     * @brief Current cost value (useful at some points e.g. linesearches)
     */
    double costCurr_;

    /**
     * Output of SSOM for R
     */
    SomUtils::VecMatD Rcurr_;

    /**
     * Output of SSOM for T
     */
    SomUtils::MatD Tcurr_;

    /**
     * Output of SSOM for Lambdas
     */
    SomUtils::MatD LambdasCurr_;

    /**
     * Output of SSOM for R
     */
    SomUtils::VecMatD Rout_;

    /**
     * Output of SSOM for T
     */
    SomUtils::MatD Tout_;

    /**
     * Output of SSOM for Lambdas
     */
    SomUtils::MatD LambdasOut_;

    /**
     * Threshold for early stopping
     */
    double transfEndThresh_;

    /**
     * Maximum number of iterations
     * One iteration includes both stepOne and stepTwo
     */
    int maxNumIterations_;
};

#endif /*SSOM_PROCRUSTES_H_*/

// int ar = A.rows();
// integer n[] = {ar};
// int br = b.rows();
// integer nrhs[] = {br};
// int ac = A.cols();
// integer lda[] = {ac};
// int bc = b.cols();
// integer ldb[] = {bc};
// int *ipiv, *info;
// dgesv_(n, nrhs, A.data(), lda, ipiv, b.data(), ldb, info);
// /* Check for the exact singularity */
// if (info[0] > 0)
// {
//     printf("The diagonal element of the triangular factor of A,\n");
//     printf("U(%i,%i) is zero, so that A is singular;\n", info[0], info[0]);
//     printf("the solution could not be computed.\n");
//     exit(1);
// }
