#ifndef SOM_PROCRUSTES_H_
#define SOM_PROCRUSTES_H_

#include "som_utils.h"
#include "thirdparty/roptlib/cwrapper/lapack/dgesv.h"

class SomProcrustes
{
public:
    /**
     * Default (empty) constructor
     */
    SomProcrustes();

    /**
     * Standard constructor with problem data inputs (most commonly used)
     */
    SomProcrustes(SomUtils::SomSize somSz, SomUtils::MatD &Tijs, Eigen::MatrixXi &edges);

    /**
     * Standard constructor with problem data inputs (most commonly used)
     */
    SomProcrustes(const SomUtils::SomSize somSz, const SomUtils::MatD &Tijs, const Eigen::MatrixXi &edges);

    /**
     * Default (empty) destructor
     */
    virtual ~SomProcrustes();

private:
    /**
     * Extract i-th part from Mmat or Nmat (see makeMN)
     */
    void getMiFromMmat(const SomUtils::MatD &Mmat, SomUtils::MatD &Mi, int idx) const;

    /**
     * Estimate rotations
     */
    void stepOne(const SomUtils::MatD &T);

    /**
     * Estimate translations
     */
    void stepTwo(const SomUtils::VecMatD &R);

    /**
     * Make matrices M and N (used in stepOne)
     */
    void makeMN(const SomUtils::MatD &T, SomUtils::MatD &Mmat, SomUtils::MatD &Nmat) const;

public:
    /**
     * Make A and b matrices (used in stepTwo)
     */
    void makeAb(const SomUtils::VecMatD &Rgf, SomUtils::MatD &A, SomUtils::MatD &b) const;

    /**
     * Run the entire pipeline
     */
    void run();

    /**
     * Set Tstart_ before running stepOne for the first time
     */
    void setTstart(const SomUtils::MatD &Tstart);

    /**
     * Set ground truth information object for R
     */
    void setRgt(const SomUtils::VecMatD &Rgt);

    /**
     * Set ground truth information object for T
     */
    void setTgt(const SomUtils::MatD &Tgt);

    /**
     * Get the output of SOM-Procrustes for T
     */
    SomUtils::MatD getTout() const;

    /**
     * Get the output of SOM-Procrustes for T
     */
    SomUtils::VecMatD getRout() const;

private:
    /**
     * Vectorize R and T into a single vector XvecOut
     * equivalent of Matla's XvecOut = [R(:); T(:)];
     */
    void vectorizeRT(const SomUtils::VecMatD &R, const SomUtils::MatD &T, SomUtils::MatD &XvecOut) const;

    /**
     * Compute norm of [R1(:); T1(:)] - [R2(:); T2(:)]
     */
    double norm(const SomUtils::VecMatD &R1, const SomUtils::MatD &T1, const SomUtils::VecMatD &R2, const SomUtils::MatD &T2) const;

    /**
     * Compute norm of [R(:); T(:)]
     */
    double norm(const SomUtils::VecMatD &R, const SomUtils::MatD &T) const;

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
    SomUtils::VecMatD Rcurr_;

    /**
     * Output of RSOM for T
     */
    SomUtils::MatD Tcurr_;

    /**
     * Output of RSOM for R
     */
    SomUtils::VecMatD Rout_;

    /**
     * Output of RSOM for T
     */
    SomUtils::MatD Tout_;

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

#endif /*SOM_PROCRUSTES_H_*/

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
