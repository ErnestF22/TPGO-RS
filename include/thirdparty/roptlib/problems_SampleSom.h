// #include "manifolds_Rotations.h"
#include "manifolds_Stiefel.h"
#include "manifolds_Euclidean.h"
#include "manifolds_MultiManifolds.h"

#include <vector>
#include <queue>
#include <numeric>

#include <eigen3/Eigen/Dense>
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
         * Compute cost function, given input x
         * OBS. ROPT to Eig conversion performed here internally
         */
        virtual realdp f(const Variable &x) const;

        double costEigenSEdN(const SomUtils::MatD &xEigen) const
        {
            double cost = 0.0f;
            for (int e = 0; e < numEdges_; ++e)
            {
                SomUtils::MatD Ri(SomUtils::MatD::Zero(sz_.d_, sz_.d_));
                SomUtils::MatD Ti(SomUtils::MatD::Zero(sz_.d_, 1));
                SomUtils::MatD Tj(SomUtils::MatD::Zero(sz_.d_, 1));

                SomUtils::VecD tij(SomUtils::VecD::Zero(sz_.d_));
                tij = Tijs_.col(e);

                int i = edges_(e, 0) - 1; // !! -1
                int j = edges_(e, 1) - 1; // !! -1
                getRiSEdN(xEigen, Ri, i);
                getRiSEdN(xEigen, Ti, i);
                getRiSEdN(xEigen, Tj, j);

                // ROFL_VAR3(i, j, e);
                // ROFL_VAR4(Ri, tij.transpose(), Ti.transpose(), Tj.transpose());

                double costE = (Ri * tij - Tj + Ti).norm(); // TODO: use squaredNorm() here directly
                // ROFL_VAR1(costE);

                cost += costE * costE;
            }
            return cost;
        }

        double costEigen(const SomUtils::MatD &xEigen) const
        {
            double cost = 0.0f;
            for (int e = 0; e < numEdges_; ++e)
            {
                SomUtils::MatD Ri(SomUtils::MatD::Zero(sz_.p_, sz_.d_));
                SomUtils::MatD Ti(SomUtils::MatD::Zero(sz_.p_, 1));
                SomUtils::MatD Tj(SomUtils::MatD::Zero(sz_.p_, 1));

                SomUtils::VecD tij(SomUtils::VecD::Zero(sz_.d_));
                tij = Tijs_.col(e);

                int i = edges_(e, 0) - 1; // !! -1
                int j = edges_(e, 1) - 1; // !! -1
                getRi(xEigen, Ri, i);
                getTi(xEigen, Ti, i);
                getTi(xEigen, Tj, j);

                // ROFL_VAR3(i, j, e);
                // ROFL_VAR4(Ri, tij.transpose(), Ti.transpose(), Tj.transpose());

                double costE = (Ri * tij - Tj + Ti).norm(); // TODO: use squaredNorm() here directly
                // ROFL_VAR1(costE);

                cost += costE * costE;
            }
            return cost;
        }

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
         * Vertically stack input vector
         */
        void vstack(const SomUtils::VecMatD &in, SomUtils::MatD &out) const;

        /**
         * Horizontally stack input vector
         */
        void hstack(const SomUtils::VecMatD &in, SomUtils::MatD &out) const;

        /**
         * Unstack a vertically stacked "3D" array
         */
        void unStackV(const SomUtils::MatD &in, SomUtils::VecMatD &out, int rowsOut = 3) const;

        /**
         * Unstack a horizontally stacked "3D" array
         */
        void unStackH(const SomUtils::MatD &in, SomUtils::VecMatD &out, int colsOut = 3) const;

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
         * Project Hin onto tangent space at Y (3D Stiefel)
         */
        void stiefelTangentProj(const SomUtils::VecMatD &Y, const SomUtils::VecMatD &Hin, SomUtils::VecMatD &Hout) const;

        /**
         * Project Hin onto tangent space at Y (single Stiefel matrix)
         */
        void stiefelTangentProj(const SomUtils::MatD &Y, const SomUtils::MatD &Hin, SomUtils::MatD &Hout) const;

        /**
         * Extract symmetric part of input (square) matrix
         */
        void extractSymmetricPart(const SomUtils::MatD &in, SomUtils::MatD &out) const;

        // private: //TODO: separate public from private members

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

        bool eigencheckHessianGenproc(const double &lambda,
                                      const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR,
                                      const SomUtils::MatD &xT, const SomUtils::MatD &uT, double thr = 1e-3) const
        {
            // if ~exist('thr','var')
            //     thr = 1e-3;
            // end

            // hessV = [matStackH(hess_fun_han(v).R), hess_fun_han(v).T];

            // OBS. TODO: These vectorizations can probably be done faster as elements order is not really important
            SomUtils::VecMatD rhR(sz_.n_, SomUtils::MatD::Zero(xR[0].rows(), xR[0].cols()));
            SomUtils::MatD rhT(SomUtils::MatD::Zero(xT.rows(), xT.cols()));
            hessGenprocEigen(xR, uR, xT, uT, rhR, rhT);

            int fullRotsSz = sz_.n_ * sz_.p_ * sz_.d_; // TODO: class function for this
            SomUtils::MatD rhRvec(SomUtils::MatD::Zero(fullRotsSz, 1));
            int fullIdx = 0;
            for (int i = 0; i < sz_.n_; ++i)
            {
                for (int j = 0; j < sz_.d_; ++j)
                {
                    for (int k = 0; k < sz_.p_; ++k)
                    {
                        rhRvec(fullIdx, 0) = rhR[i](k, j);
                        fullIdx++;
                    }
                }
            }
            int fullTranslSz = sz_.n_ * sz_.p_; // TODO: class member function for this
            SomUtils::MatD rhTvec(SomUtils::MatD::Zero(fullTranslSz, 1));
            fullIdx = 0; //!! resetting fullIdx
            for (int j = 0; j < sz_.n_; ++j)
            {
                for (int k = 0; k < sz_.p_; ++k)
                {
                    rhTvec(fullIdx, 0) = rhT(k, j);
                    fullIdx++;
                }
            }
            Eigen::VectorXd hessV(rhRvec.size() + rhTvec.size());
            hessV << rhRvec, rhTvec;

            // uFull = [matStackH(v.R), v.T];
            SomUtils::MatD uRvec(SomUtils::MatD::Zero(fullRotsSz, 1));
            fullIdx = 0;
            for (int i = 0; i < sz_.n_; ++i)
            {
                for (int j = 0; j < sz_.d_; ++j)
                {
                    for (int k = 0; k < sz_.p_; ++k)
                    {
                        uRvec(fullIdx, 0) = uR[i](k, j);
                        fullIdx++;
                    }
                }
            }
            SomUtils::MatD uTvec(SomUtils::MatD::Zero(fullTranslSz, 1));
            fullIdx = 0; //!! resetting fullIdx
            for (int j = 0; j < sz_.n_; ++j)
            {
                for (int k = 0; k < sz_.p_; ++k)
                {
                    uTvec(fullIdx, 0) = uT(k, j);
                    fullIdx++;
                }
            }
            Eigen::VectorXd uFull(uRvec.size() + uTvec.size());
            uFull << uRvec, uTvec;

            // diff = norm((lambda)*uFull(:) - hessV(:),'inf');
            double diff = (lambda * uFull - hessV).lpNorm<Eigen::Infinity>();

            // if diff < thr
            //     fprintf("%g is a GENPROC eigenvalue\n", lambda);
            //     eig_bool = boolean(1);
            // else
            //     fprintf("%g is NOT a GENPROC eigenvalue: diff %g > thr\n", ...
            //         lambda, diff);
            //     eig_bool = boolean(0);
            // end

            if (diff < thr)
            {
                std::cout << lambda << " is a GENPROC eigenvalue" << std::endl;
                return true;
            }
            else
            {
                std::cout << lambda << " is NOT a GENPROC eigenvalue: diff " << diff << " > thr " << thr << std::endl;
                // return false;
            }

            return false;
        }

        bool eigencheckHessianGenprocShifted(const double &lambda,
                                             const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR,
                                             const SomUtils::MatD &xT, const SomUtils::MatD &uT, double mu, double thr = 1e-3) const
        {
            // if ~exist('thr','var')
            //     thr = 1e-3;
            // end

            // hessV = [matStackH(hess_fun_han(v).R), hess_fun_han(v).T];

            // OBS. TODO: These vectorizations can probably be done faster as elements order is not really important
            SomUtils::VecMatD rhR(sz_.n_, SomUtils::MatD::Zero(xR[0].rows(), xR[0].cols()));
            SomUtils::MatD rhT(SomUtils::MatD::Zero(xT.rows(), xT.cols()));
            hessGenprocEigenShifted(xR, uR, xT, uT, mu, rhR, rhT);

            int fullRotsSz = sz_.n_ * sz_.p_ * sz_.d_; // TODO: class function for this
            SomUtils::MatD rhRvec(SomUtils::MatD::Zero(fullRotsSz, 1));
            int fullIdx = 0;
            for (int i = 0; i < sz_.n_; ++i)
            {
                for (int j = 0; j < sz_.d_; ++j)
                {
                    for (int k = 0; k < sz_.p_; ++k)
                    {
                        rhRvec(fullIdx, 0) = rhR[i](k, j);
                        fullIdx++;
                    }
                }
            }
            int fullTranslSz = sz_.n_ * sz_.p_; // TODO: class member function for this
            SomUtils::MatD rhTvec(SomUtils::MatD::Zero(fullTranslSz, 1));
            fullIdx = 0; //!! resetting fullIdx
            for (int j = 0; j < sz_.n_; ++j)
            {
                for (int k = 0; k < sz_.p_; ++k)
                {
                    rhTvec(fullIdx, 0) = rhT(k, j);
                    fullIdx++;
                }
            }
            Eigen::VectorXd hessV(rhRvec.size() + rhTvec.size());
            hessV << rhRvec, rhTvec;

            // uFull = [matStackH(v.R), v.T];
            SomUtils::MatD uRvec(SomUtils::MatD::Zero(fullRotsSz, 1));
            fullIdx = 0;
            for (int i = 0; i < sz_.n_; ++i)
            {
                for (int j = 0; j < sz_.d_; ++j)
                {
                    for (int k = 0; k < sz_.p_; ++k)
                    {
                        uRvec(fullIdx, 0) = uR[i](k, j);
                        fullIdx++;
                    }
                }
            }
            SomUtils::MatD uTvec(SomUtils::MatD::Zero(fullTranslSz, 1));
            fullIdx = 0; //!! resetting fullIdx
            for (int j = 0; j < sz_.n_; ++j)
            {
                for (int k = 0; k < sz_.p_; ++k)
                {
                    uTvec(fullIdx, 0) = uT(k, j);
                    fullIdx++;
                }
            }
            Eigen::VectorXd uFull(uRvec.size() + uTvec.size());
            uFull << uRvec, uTvec;

            // diff = norm((lambda)*uFull(:) - hessV(:),'inf');
            double diff = (lambda * uFull - hessV).lpNorm<Eigen::Infinity>();

            // if diff < thr
            //     fprintf("%g is a GENPROC eigenvalue\n", lambda);
            //     eig_bool = boolean(1);
            // else
            //     fprintf("%g is NOT a GENPROC eigenvalue: diff %g > thr\n", ...
            //         lambda, diff);
            //     eig_bool = boolean(0);
            // end

            if (diff < thr)
            {
                std::cout << lambda << " is a GENPROC eigenvalue" << std::endl;
                return true;
            }
            else
            {
                std::cout << lambda << " is NOT a GENPROC eigenvalue: diff " << diff << " > thr " << thr << std::endl;
                // return false;
            }

            return false;
        }

        void catZeroRow(const SomUtils::MatD &mIn, SomUtils::MatD &mOut) const
        {
            ROFL_ASSERT(mOut.rows() == mIn.rows() + 1);
            ROFL_ASSERT(mOut.cols() == mIn.cols());

            mOut.setZero();
            mOut.block(0, 0, mIn.rows(), mIn.cols()) = mIn;
        }

        void catZeroRow3dArray(const SomUtils::VecMatD &mIn, SomUtils::VecMatD &mOut) const
        {
            ROFL_ASSERT(mIn.size() == mOut.size())

            int n = mIn.size();
            for (int i = 0; i < n; ++i)
            {
                catZeroRow(mIn[i], mOut[i]);
            }
        }

        void normalizeEucl(const SomUtils::MatD &mIn, SomUtils::MatD &mOut) const
        {
            ROFL_ASSERT(mIn.rows() == mOut.rows() && mIn.cols() == mOut.cols())

            mOut = mIn;
            double normF = mIn.norm(); // TODO: maybe use .normalized() directly?

            mOut /= normF;
            // ROFL_VAR3(mIn, normF, mOut);
        }

        void normalizeEucl(const SomUtils::VecMatD &mIn, SomUtils::VecMatD &mOut) const
        {
            ROFL_ASSERT(mIn.size() == mOut.size())
            std::for_each(mOut.begin(), mOut.end(), [](SomUtils::MatD &x) { //^^^ take argument by reference: LAMBDA FUNCTION
                x.setZero();
            });

            int n = mIn.size();
            for (int i = 0; i < n; ++i)
                normalizeEucl(mIn[i], mOut[i]);
        }

        void RoptToEigStiefel(Vector x, SomUtils::MatD &xEigen) const
        {
            Vector xT = x.GetTranspose(); // Eigen ADV init is row-major!!

            int totSz = xEigen.rows(); // xEigen is supposed to be a vectorized matrix

            const realdp *xArr = xT.ObtainWriteEntireData();
            for (int i = 0; i < totSz; ++i)
                xEigen(i) = xArr[i];
        }

        void stiefelRandTgNormVector(const SomUtils::MatD &mIn, SomUtils::MatD &mOut) const
        {
            mOut.setIdentity();

            int r = mIn.rows();
            int c = mIn.cols();

            ROFL_ASSERT(r == mOut.rows() && c == mOut.cols());

            // generate random Stiefel element
            Vector tmp = Stiefel(r, c).RandominManifold();
            // tmp.Print("tmp inside stiefelRandTgNormVector()");
            SomUtils::MatD tmpEigVec(SomUtils::MatD::Zero(r * c, 1));
            RoptToEigStiefel(tmp, tmpEigVec); // xEigen is supposed to be a vectorized matrix
            // ROFL_VAR1(tmpEigVec);
            SomUtils::MatD tmpEig = tmpEigVec.reshaped<Eigen::RowMajor>(r, c);
            // ROFL_VAR1(tmpEig);
            SomUtils::MatD tmpProj(SomUtils::MatD::Zero(r, c));
            stiefelTangentProj(mIn, tmpEig, tmpProj);
            // ROFL_VAR1(tmpProj);

            normalizeEucl(tmpProj, mOut);
            // ROFL_VAR1(mOut);

            // For any matrix representative $U \in St(n, p)$, the tangent space of $St(n, p)$ at $U$ is represented by
            // U\transpose \Delta = -\Delta\transpose U

            double diff = (mIn.transpose() * mOut + mOut.transpose() * mIn).cwiseAbs().maxCoeff();
            ROFL_ASSERT_VAR1(diff >= 0 && diff < 1e-5, diff);
            double mOutNorm = mOut.norm();
            ROFL_ASSERT_VAR1(mOutNorm > 1 - 1e-5 && mOutNorm < 1 + 1e-5, mOutNorm);
        }

        void stiefelRandTgNormVector(const SomUtils::VecMatD &mIn, SomUtils::VecMatD &mOut) const
        {
            int n = mIn.size();
            ROFL_ASSERT(n == mOut.size());
            std::for_each(mOut.begin(), mOut.end(), [](SomUtils::MatD &x) { //^^^ take argument by reference: LAMBDA FUNCTION
                x.setZero();
            });

            for (int i = 0; i < n; ++i)
            {
                stiefelRandTgNormVector(mIn[i], mOut[i]);
            }
        }

        void vectorizeR(const SomUtils::VecMatD &R, SomUtils::MatD &RvecOut) const
        {
            // int fullRotsSz = sz_.p_ * sz_.d_ * sz_.n_;

            // for (int i=0; i<fullRotsSz; ++i) {
            // }

            int n = R.size();

            int fullIdx = 0;
            for (int i = 0; i < n; ++i)
            {
                int ric = R[i].cols();
                int rir = R[i].rows();
                for (int j = 0; j < ric; ++j)
                {
                    for (int k = 0; k < rir; ++k)
                    {
                        RvecOut(fullIdx, 0) = R[i](k, j);
                        fullIdx++;
                        // ROFL_VAR4(i, j, k, fullIdx);
                    }
                }
            }
            ROFL_ASSERT(fullIdx == RvecOut.rows())
        }

        void pimFunctionGenproc(const SomUtils::VecMatD &xR, const SomUtils::MatD &xT, const SomUtils::VecMatD &uR, const SomUtils::MatD &uT,
                                double &lambdaMax, SomUtils::VecMatD &uOutR, SomUtils::MatD &uOutT, double thresh = 1e-5) const
        {
            // Note: normalization is done across entire ProdMani vector through simple eucl. metric

            // % % R iterative_change = 1e+6;
            // xR = x_start.R;
            // xT = x_start.T;
            // xfull = [ matStackH(x_start.R), x_start.T ];

            int staircaseLevel = xT.rows();
            ROFL_VAR1(staircaseLevel);

            SomUtils::MatD uRhStacked(SomUtils::MatD::Zero(staircaseLevel, sz_.d_ * sz_.n_));
            // ROFL_VAR1("hstack call from here");
            hstack(uR, uRhStacked);
            ROFL_VAR1(uRhStacked);

            SomUtils::MatD uTcopy = uT; // useful for keeping const in function params
            ROFL_VAR1(uTcopy);

            SomUtils::MatD uFullHst(SomUtils::MatD::Zero(staircaseLevel, sz_.n_ + uRhStacked.cols()));
            uFullHst.block(0, 0, staircaseLevel, uRhStacked.cols()) = uRhStacked;
            uFullHst.block(0, uRhStacked.cols(), staircaseLevel, sz_.n_) = uTcopy;
            ROFL_VAR1(uFullHst);

            // iteration_num = 0;
            int iterationNum = 0;
            // while (iteration_num < 2500)
            //     % &&(abs(iterative_change) > thresh)
            //      iteration_num = iteration_num + 1;
            //      x_prev_R = xR;
            //      x_prev_T = xT;
            //      xfull_prev = [ matStackH(x_prev_R), x_prev_T ];
            //      norm_RT = norm(xfull);
            //      x.R = xR / norm_RT;
            //      x.T = xT / norm_RT;
            //      fx = f(x);
            //      xR = -fx.R;
            //      xT = -fx.T;
            //      xfull = [ matStackH(xR), xT ];
            //      iterative_change = max(normalization_fun(xfull_prev - xfull), [], "all");
            // end
            SomUtils::MatD uRprevHst(SomUtils::MatD::Zero(staircaseLevel, sz_.d_ * sz_.n_));
            SomUtils::MatD uTprev(SomUtils::MatD::Zero(staircaseLevel, sz_.n_));
            while (iterationNum < 2500) // && iterativeChange < 1e-3
            {
                // ROFL_VAR1(iterationNum);
                iterationNum++;

                // 1
                uRprevHst = uRhStacked;
                // ROFL_VAR1(uRprevHst);

                uTprev = uTcopy;
                // ROFL_VAR1(uTprev);

                // 2
                SomUtils::MatD uFullHstPrev(SomUtils::MatD::Zero(staircaseLevel, sz_.n_ + uRhStacked.cols()));
                uFullHstPrev.block(0, 0, staircaseLevel, uRprevHst.cols()) = uRprevHst;
                uFullHstPrev.block(0, uRprevHst.cols(), staircaseLevel, uTprev.cols()) = uTprev;
                // ROFL_VAR1(uFullHstPrev);

                // 3
                double normRT = uFullHst.norm();
                // ROFL_VAR1(normRT);

                // 4
                uRhStacked /= normRT;
                uTcopy /= normRT;

                SomUtils::VecMatD uRunstackedTmp(sz_.n_, SomUtils::MatD::Zero(staircaseLevel, sz_.d_));
                unStackH(uRhStacked, uRunstackedTmp);

                // 5
                SomUtils::VecMatD uRunstackedOutTmp(sz_.n_, SomUtils::MatD::Zero(staircaseLevel, sz_.d_));
                SomUtils::MatD uTout(SomUtils::MatD::Zero(staircaseLevel, sz_.n_));
                hessGenprocEigen(xR, uRunstackedTmp, xT, uTcopy, uRunstackedOutTmp, uTout);

                // 6
                std::for_each(uRunstackedOutTmp.begin(), uRunstackedOutTmp.end(), [](SomUtils::MatD &x) { //^^^ take argument by reference: LAMBDA FUNCTION
                    x *= -1;
                });
                uTcopy = -uTout;

                // 7
                hstack(uRunstackedOutTmp, uRhStacked);
                // ROFL_VAR1("hstack call from here");
                uFullHst.block(0, 0, staircaseLevel, uRhStacked.cols()) = uRhStacked;
                uFullHst.block(0, uRhStacked.cols(), staircaseLevel, sz_.n_) = uTcopy;
                // ROFL_VAR1(uFullHst);

                //      iterative_change = max(normalization_fun(xfull_prev - xfull), [], "all");
                double iterativeChange = (uFullHstPrev - uFullHst).cwiseAbs().maxCoeff();
                ROFL_VAR2(iterationNum, iterativeChange);
            }

            // 1
            // norm_RT_max = norm([ matStackH(x.R), x.T ]);
            double normRTmax = uFullHst.norm();

            // 2
            // x_max.R = x.R / norm_RT_max;
            // x_max.T = x.T / norm_RT_max;
            uFullHst /= normRTmax;

            // 3
            // f_x_max = f(x_max);
            SomUtils::VecMatD uRunstacked(sz_.n_, SomUtils::MatD::Zero(staircaseLevel, sz_.d_));
            unStackH(uRhStacked / normRTmax, uRunstacked, sz_.d_);
            SomUtils::VecMatD fxRunstackedOut(sz_.n_, SomUtils::MatD::Zero(staircaseLevel, sz_.d_));

            SomUtils::MatD uTout1(SomUtils::MatD::Zero(staircaseLevel, sz_.n_));

            hessGenprocEigen(xR, uRunstacked, xT, uTcopy / normRTmax, fxRunstackedOut, uTout1);

            // 4
            SomUtils::MatD vecR(SomUtils::MatD::Zero(sz_.d_ * staircaseLevel * sz_.n_, 1));
            vectorizeR(uRunstacked, vecR);
            SomUtils::MatD vecFxR(SomUtils::MatD::Zero(sz_.d_ * staircaseLevel * sz_.n_, 1));
            vectorizeR(fxRunstackedOut, vecFxR);
            // lambda_max = x_max.R( :) ' * f_x_max.R(:) + x_max.T(:)' * f_x_max.T( :);
            auto lmax = vecR.transpose() * vecFxR + uTout1.reshaped(1, staircaseLevel * sz_.n_) * (uTcopy / normRTmax).reshaped(staircaseLevel * sz_.n_, 1);
            ROFL_VAR1(lmax)
            lambdaMax = lmax(0, 0); // lmax is supposedly a scalar

            // 5
            // full_xmax = [ matStackH(x_max.R), x_max.T ];
            auto uFullRhSt = uFullHst.block(0, 0, staircaseLevel, sz_.d_ * sz_.n_);
            unStackH(uFullRhSt, uOutR, sz_.d_);
            uOutT = uFullHst.block(0, sz_.d_ * sz_.n_, staircaseLevel, sz_.n_);

            // 6
            // lambda_max = lambda_max / sum(stiefel_metric([], full_xmax( :), full_xmax( :)));
            lambdaMax /= uFullHst.norm();
        }

        void pimFunctionGenprocShifted(const SomUtils::VecMatD &xR, const SomUtils::MatD &xT, const SomUtils::VecMatD &uR, const SomUtils::MatD &uT, double mu,
                                       double &lambdaMax, SomUtils::VecMatD &uOutR, SomUtils::MatD &uOutT, double thresh = 1e-5) const
        {
            // Note: normalization is done across entire ProdMani vector through simple eucl. metric

            // % % R iterative_change = 1e+6;
            // xR = x_start.R;
            // xT = x_start.T;
            // xfull = [ matStackH(x_start.R), x_start.T ];

            int staircaseLevel = xT.rows();

            SomUtils::MatD uRhStacked(SomUtils::MatD::Zero(staircaseLevel, sz_.d_ * sz_.n_));
            // ROFL_VAR1("hstack call from here");
            hstack(uR, uRhStacked);

            SomUtils::MatD uTcopy = uT; // useful for keeping const in function params

            SomUtils::MatD uFullHst(SomUtils::MatD::Zero(staircaseLevel, sz_.n_ + uRhStacked.cols()));
            uFullHst.block(0, 0, staircaseLevel, uRhStacked.cols()) = uRhStacked;
            uFullHst.block(0, uRhStacked.cols(), staircaseLevel, sz_.n_) = uTcopy;
            // iteration_num = 0;
            int iterationNum = 0;
            // while (iteration_num < 2500)
            //     % &&(abs(iterative_change) > thresh)
            //      iteration_num = iteration_num + 1;
            //      x_prev_R = xR;
            //      x_prev_T = xT;
            //      xfull_prev = [ matStackH(x_prev_R), x_prev_T ];
            //      norm_RT = norm(xfull);
            //      x.R = xR / norm_RT;
            //      x.T = xT / norm_RT;
            //      fx = f(x);
            //      xR = -fx.R;
            //      xT = -fx.T;
            //      xfull = [ matStackH(xR), xT ];
            //      iterative_change = max(normalization_fun(xfull_prev - xfull), [], "all");
            // end
            SomUtils::MatD uRprevHst(SomUtils::MatD::Zero(staircaseLevel, sz_.d_ * sz_.n_));
            SomUtils::MatD uTprev(SomUtils::MatD::Zero(staircaseLevel, sz_.n_));
            while (iterationNum < 2500) // && iterativeChange < 1e-3
            {
                iterationNum++;
                uRprevHst = uRhStacked;
                uTprev = uTcopy;

                SomUtils::MatD uFullHstPrev(SomUtils::MatD::Zero(staircaseLevel, sz_.n_ + uRhStacked.cols()));
                uFullHstPrev.block(0, 0, staircaseLevel, uRhStacked.cols()) = uRprevHst;
                uFullHstPrev.block(0, uRhStacked.cols(), staircaseLevel, sz_.n_) = uTprev;

                double normRT = uFullHst.norm();
                uRhStacked /= normRT;
                uTcopy /= normRT;

                SomUtils::VecMatD uRunstackedTmp(sz_.n_, SomUtils::MatD::Zero(staircaseLevel, sz_.d_));
                unStackH(uRhStacked, uRunstackedTmp);

                SomUtils::VecMatD uRunstackedOutTmp(sz_.n_, SomUtils::MatD::Zero(staircaseLevel, sz_.d_));
                SomUtils::MatD uTout(SomUtils::MatD::Zero(staircaseLevel, sz_.n_));
                hessGenprocEigenShifted(xR, uRunstackedTmp, xT, uTcopy, mu, uRunstackedOutTmp, uTout);

                std::for_each(uRunstackedOutTmp.begin(), uRunstackedOutTmp.end(), [](SomUtils::MatD &x) { //^^^ take argument by reference: LAMBDA FUNCTION
                    x *= -1;
                });
                uTcopy = -uTout;

                // ROFL_VAR1("hstack call from here");
                hstack(uRunstackedOutTmp, uRhStacked);
                uFullHst.block(0, 0, staircaseLevel, uRhStacked.cols()) = uRhStacked;
                uFullHst.block(0, uRhStacked.cols(), staircaseLevel, sz_.n_) = uTcopy;

                //      iterative_change = max(normalization_fun(xfull_prev - xfull), [], "all");
                double iterativeChange = (uFullHstPrev - uFullHst).cwiseAbs().maxCoeff();
                ROFL_VAR2(iterationNum, iterativeChange);
            }

            // norm_RT_max = norm([ matStackH(x.R), x.T ]);
            double normRTmax = uFullHst.norm();

            // x_max.R = x.R / norm_RT_max;
            // x_max.T = x.T / norm_RT_max;
            uFullHst /= normRTmax;

            // f_x_max = f(x_max);
            SomUtils::VecMatD uRunstacked(sz_.n_, SomUtils::MatD::Zero(staircaseLevel, sz_.d_));
            unStackH(uRhStacked / normRTmax, uRunstacked, sz_.d_);
            SomUtils::VecMatD fxRunstackedOut(sz_.n_, SomUtils::MatD::Zero(staircaseLevel, sz_.d_));

            SomUtils::MatD uTout1(SomUtils::MatD::Zero(staircaseLevel, sz_.n_));

            hessGenprocEigenShifted(xR, uRunstacked, xT, uTcopy / normRTmax, mu, fxRunstackedOut, uTout1);

            // % lambda_max_R = sum(stiefel_metric([], (x_max.R), f_x_max.R)) / ... % sum(stiefel_metric([], x_max.R, x_max.R));
            // % lambda_max_T = sum(stiefel_metric([], (x_max.T), f_x_max.T)) / ... % sum(stiefel_metric([], x_max.T, x_max.T));

            SomUtils::MatD vecR(SomUtils::MatD::Zero(sz_.d_ * staircaseLevel * sz_.n_, 1));
            vectorizeR(uRunstacked, vecR);
            SomUtils::MatD vecFxR(SomUtils::MatD::Zero(sz_.d_ * staircaseLevel * sz_.n_, 1));
            vectorizeR(fxRunstackedOut, vecFxR);

            // lambda_max = x_max.R( :) ' * f_x_max.R(:) + x_max.T(:)' * f_x_max.T( :);
            auto lmax = vecR.transpose() * vecFxR + uTout1.reshaped(1, staircaseLevel * sz_.n_) * (uTcopy / normRTmax).reshaped(staircaseLevel * sz_.n_, 1);
            lambdaMax = lmax(0, 0);

            // full_xmax = [ matStackH(x_max.R), x_max.T ];
            auto uFullRhSt = uFullHst.block(0, 0, staircaseLevel, sz_.d_ * sz_.n_);
            unStackH(uFullRhSt, uOutR, sz_.d_);
            uOutT = uFullHst.block(0, sz_.d_ * sz_.n_, staircaseLevel, sz_.n_);

            // lambda_max = lambda_max / sum(stiefel_metric([], full_xmax( :), full_xmax( :)));
            lambdaMax /= uFullHst.norm();
        }

        void rsomPimHessianGenprocSmall(double thresh, const SomUtils::VecMatD &R, const SomUtils::MatD &T) const
        {
            // [Y_star, lambda, v] = rsom_pim_hessian_genproc( ...
            //     X, problem_struct_next, thr);
            // disp("v") // %just to remove unused variable warning
            // disp(v)
            // if lambda > 0
            //     disp("R, T eigenvals > 0: exiting staircase")
            // break;

            /////////////////////////////////////////////////////
            // if ~exist('thresh', 'var')
            //     thresh = 1e-6;
            // end

            // Rnext = cat_zero_rows_3d_array(X.R);
            // Tnext = cat_zero_row(X.T);
            // Xnext.R = Rnext;
            // Xnext.T = Tnext;
            // rhess_fun_han = @(u) hess_genproc(Xnext,u,problem_struct_next);
            int staircaseNextStepLevel = T.rows() + 1;
            ROFL_VAR1(staircaseNextStepLevel);
            SomUtils::VecMatD Rnext(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            catZeroRow3dArray(R, Rnext);
            SomUtils::MatD Tnext(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
            catZeroRow(T, Tnext);

            // stiefel_normalize_han = @(x) x./ (norm(x(:))); //Note: this is basically eucl_normalize_han

            // u_start.R = stiefel_randTangentNormVector(Rnext);
            SomUtils::VecMatD RnextTg(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            stiefelRandTgNormVector(Rnext, RnextTg);
            // u_start.R = stiefel_normalize(Rnext, u_start.R);
            SomUtils::VecMatD RnextTgNorm(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            normalizeEucl(RnextTg, RnextTgNorm);

            // u_start.T = rand(size(Tnext));
            auto TnextTg = SomUtils::MatD::Random(staircaseNextStepLevel, sz_.n_);
            // u_start.T = stiefel_normalize_han(u_start.T);
            SomUtils::MatD TnextTgNorm(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
            normalizeEucl(TnextTg, TnextTgNorm);

            // [lambda_pim, v_pim] = pim_function_genproc(rhess_fun_han, u_start, stiefel_normalize_han, thresh);
            // disp('Difference between lambda*v_max and H(v_max) should be in the order of the tolerance:')
            double lambdaPim = 1e+6;
            SomUtils::VecMatD vPimR(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            SomUtils::MatD vPimT(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));

            pimFunctionGenproc(Rnext, Tnext, RnextTgNorm, TnextTgNorm, lambdaPim, vPimR, vPimT);
            std::cout << "Difference between lambda_pim_after_shift*v_pim_after_shift"
                      << " and H_SH(v_pim_after_shift) should be in the order of the tolerance:" << std::endl;
            eigencheckHessianGenproc(lambdaPim, Rnext, vPimR, Tnext, vPimT);

            // if lambda_pim>0
            double highestNormEigenval = 1e+6;
            if (lambdaPim > 0)
            {
                std::cout << "lambdaPim " << lambdaPim << std::endl;
                double mu = 1.1 * lambdaPim;

                //     rhess_shifted_fun_han = @(u) hess_genproc_shifted(Xnext,u,mu,problem_struct_next);

                //     // %run shifted power iteration
                //     u_start_second_iter.R = stiefel_randTangentNormVector(Rnext);
                //     u_start_second_iter.R = stiefel_normalize(Rnext, u_start_second_iter.R);
                //     u_start_second_iter.T = rand(size(Tnext));
                //     u_start_second_iter.T = stiefel_normalize_han(u_start.T);
                //     [lambda_pim_after_shift, v_pim_after_shift] = pim_function_genproc( ...
                //         rhess_shifted_fun_han, u_start_second_iter, stiefel_normalize_han, thresh);
                SomUtils::VecMatD RnextTgShift(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
                stiefelRandTgNormVector(Rnext, RnextTgShift);
                SomUtils::VecMatD RnextTgNormShift(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
                normalizeEucl(RnextTgShift, RnextTgNormShift);

                auto TnextTgShift = SomUtils::MatD::Random(staircaseNextStepLevel, sz_.n_);
                SomUtils::MatD TnextTgNormShift(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
                normalizeEucl(TnextTgShift, TnextTgNormShift);

                SomUtils::VecMatD vPimRshift(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
                SomUtils::MatD vPimTshift(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
                double lambdaPimShift = 1e+6; // "after" shift is intended
                pimFunctionGenprocShifted(Rnext, Tnext, RnextTgNorm, TnextTgNorm, mu, lambdaPimShift, vPimRshift, vPimTshift);

                //     disp(['Difference between lambda_pim_after_shift*v_pim_after_shift ' ...
                //         'and H_SH(v_pim_after_shift) should be in the order of the tolerance:'])
                //     eigencheck_hessian_genproc(lambda_pim_after_shift, v_pim_after_shift, ...
                //         rhess_shifted_fun_han);
                std::cout << "Difference between lambda_pim_after_shift*v_pim_after_shift"
                          << "and H_SH(v_pim_after_shift) should be in the order of the tolerance:" << std::endl;
                eigencheckHessianGenprocShifted(lambdaPimShift, Rnext, vPimRshift, Tnext, vPimTshift, mu);
                highestNormEigenval = lambdaPimShift + mu;
                std::cout << "Difference between (lambda_pim_after_shift + mu)*v_pim_after_shift"
                          << "and H_SH(v_pim_after_shift) should be in the order of the tolerance:" << std::endl;
                eigencheckHessianGenproc(highestNormEigenval, Rnext, vPimRshift, Tnext, vPimTshift, mu);
                // ROFL_VAR3(highestNormEigenval, vPimRshift[0], vPimTshift);
                vPimR = vPimRshift;
                vPimT = vPimTshift;
            }
            else
            {
                highestNormEigenval = lambdaPim;
            }
        }

        void rsomPimHessianGenproc(double thresh, const SomUtils::VecMatD &R, const SomUtils::MatD &T, Vector &Y0, bool armijo = true) const
        {
            // [Y_star, lambda, v] = rsom_pim_hessian_genproc( ...
            //     X, problem_struct_next, thr);
            // disp("v") // %just to remove unused variable warning
            // disp(v)
            // if lambda > 0
            //     disp("R, T eigenvals > 0: exiting staircase")
            // break;

            /////////////////////////////////////////////////////
            // if ~exist('thresh', 'var')
            //     thresh = 1e-6;
            // end

            // Rnext = cat_zero_rows_3d_array(X.R);
            // Tnext = cat_zero_row(X.T);
            // Xnext.R = Rnext;
            // Xnext.T = Tnext;
            // rhess_fun_han = @(u) hess_genproc(Xnext,u,problem_struct_next);
            int staircaseNextStepLevel = T.rows() + 1;
            ROFL_VAR1(staircaseNextStepLevel);
            SomUtils::VecMatD Rnext(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            catZeroRow3dArray(R, Rnext);
            SomUtils::MatD Tnext(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
            catZeroRow(T, Tnext);

            // stiefel_normalize_han = @(x) x./ (norm(x(:))); //Note: this is basically eucl_normalize_han

            // u_start.R = stiefel_randTangentNormVector(Rnext);
            SomUtils::VecMatD RnextTg(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            stiefelRandTgNormVector(Rnext, RnextTg);
            // u_start.R = stiefel_normalize(Rnext, u_start.R);
            SomUtils::VecMatD RnextTgNorm(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            normalizeEucl(RnextTg, RnextTgNorm);

            // u_start.T = rand(size(Tnext));
            auto TnextTg = SomUtils::MatD::Random(staircaseNextStepLevel, sz_.n_);
            // u_start.T = stiefel_normalize_han(u_start.T);
            SomUtils::MatD TnextTgNorm(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
            normalizeEucl(TnextTg, TnextTgNorm);

            // [lambda_pim, v_pim] = pim_function_genproc(rhess_fun_han, u_start, stiefel_normalize_han, thresh);
            // disp('Difference between lambda*v_max and H(v_max) should be in the order of the tolerance:')
            double lambdaPim = 1e+6;
            SomUtils::VecMatD vPimR(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            SomUtils::MatD vPimT(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));

            pimFunctionGenproc(Rnext, Tnext, RnextTgNorm, TnextTgNorm, lambdaPim, vPimR, vPimT);
            std::cout << "Difference between lambda_pim_after_shift*v_pim_after_shift"
                      << " and H_SH(v_pim_after_shift) should be in the order of the tolerance:" << std::endl;
            eigencheckHessianGenproc(lambdaPim, Rnext, vPimR, Tnext, vPimT);

            // if lambda_pim>0
            double highestNormEigenval = 1e+6;
            if (lambdaPim > 0)
            {
                std::cout << "lambdaPim " << lambdaPim << std::endl;
                double mu = 1.1 * lambdaPim;

                //     rhess_shifted_fun_han = @(u) hess_genproc_shifted(Xnext,u,mu,problem_struct_next);

                //     // %run shifted power iteration
                //     u_start_second_iter.R = stiefel_randTangentNormVector(Rnext);
                //     u_start_second_iter.R = stiefel_normalize(Rnext, u_start_second_iter.R);
                //     u_start_second_iter.T = rand(size(Tnext));
                //     u_start_second_iter.T = stiefel_normalize_han(u_start.T);
                //     [lambda_pim_after_shift, v_pim_after_shift] = pim_function_genproc( ...
                //         rhess_shifted_fun_han, u_start_second_iter, stiefel_normalize_han, thresh);
                SomUtils::VecMatD RnextTgShift(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
                stiefelRandTgNormVector(Rnext, RnextTgShift);
                SomUtils::VecMatD RnextTgNormShift(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
                normalizeEucl(RnextTgShift, RnextTgNormShift);

                auto TnextTgShift = SomUtils::MatD::Random(staircaseNextStepLevel, sz_.n_);
                SomUtils::MatD TnextTgNormShift(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
                normalizeEucl(TnextTgShift, TnextTgNormShift);

                SomUtils::VecMatD vPimRshift(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
                SomUtils::MatD vPimTshift(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
                double lambdaPimShift = 1e+6; // "after" shift is intended
                pimFunctionGenprocShifted(Rnext, Tnext, RnextTgNorm, TnextTgNorm, mu, lambdaPimShift, vPimRshift, vPimTshift);

                //     disp(['Difference between lambda_pim_after_shift*v_pim_after_shift ' ...
                //         'and H_SH(v_pim_after_shift) should be in the order of the tolerance:'])
                //     eigencheck_hessian_genproc(lambda_pim_after_shift, v_pim_after_shift, ...
                //         rhess_shifted_fun_han);
                std::cout << "Difference between lambda_pim_after_shift*v_pim_after_shift"
                          << "and H_SH(v_pim_after_shift) should be in the order of the tolerance:" << std::endl;
                eigencheckHessianGenprocShifted(lambdaPimShift, Rnext, vPimRshift, Tnext, vPimTshift, mu);
                highestNormEigenval = lambdaPimShift + mu;
                std::cout << "Difference between (lambda_pim_after_shift + mu)*v_pim_after_shift"
                          << "and H_SH(v_pim_after_shift) should be in the order of the tolerance:" << std::endl;
                eigencheckHessianGenproc(highestNormEigenval, Rnext, vPimRshift, Tnext, vPimTshift, mu);
                // ROFL_VAR3(highestNormEigenval, vPimRshift[0], vPimTshift);
                vPimR = vPimRshift;
                vPimT = vPimTshift;
            }
            else
            {
                highestNormEigenval = lambdaPim;
            } // SMALL VERSION UP TO HERE!!

            //////!!!!/////!!!!
            // // %Preparing linesearch
            // nrs_next = problem_struct_next.sz(1);
            // d = problem_struct_next.sz(2);
            // N = problem_struct_next.sz(3);
            SomUtils::SomSize szNext(staircaseNextStepLevel, sz_.d_, sz_.n_);

            Stiefel mani1(staircaseNextStepLevel, szNext.d_);
            mani1.ChooseParamsSet2();
            integer numoftypes = 2; // 2 i.e. (3D) Stiefel + Euclidean

            integer numofmani1 = szNext.n_; // num of Stiefel manifolds
            integer numofmani2 = 1;
            Euclidean mani2(staircaseNextStepLevel, szNext.n_);
            ProductManifold ProdManiNext(numoftypes, &mani1, numofmani1, &mani2, numofmani2);

            Vector xIn = ProdManiNext.RandominManifold(); //!! in other cases xIn would have been a pointer

            // tuple_next.R = stiefelfactory(nrs_next, d, N);
            // tuple_next.T = euclideanfactory(nrs_next, N);
            // M = productmanifold(tuple_next);
            // step2.M = M;
            // step2.sz = [nrs_next, d, N];
            // step2.cost = @(x) cost_genproc(x, problem_struct_next);
            // step2.grad = @(x) grad_genproc(x, problem_struct_next);
            // step2.hess = @(x, u) hess_genproc(x, u, problem_struct_next);

            // // % alpha = min(lambdas_moved) + lambdas_max;
            // // % alpha_linesearch = 10; // %TODO: set this correctly
            // // % SDPLRval = 10; // %TODO: set this correctly

            // disp("Now performing linesearch...");
            // // %Note: first output param of linesearch() would be "stepsize"
            // [~, Y0] = linesearch_decrease(step2, ...
            //     Xnext, v_pim_after_shift, cost_genproc(Xnext,problem_struct_next));

            // lambda_pim_out = highest_norm_eigenval;
            // v_pim_out = v_pim_after_shift;

            if (armijo)
            {
                { // EigToRopt scope for xIn

                    int rotSz = szNext.p_ * szNext.d_;
                    // int translSz = szNext.p_;

                    int gElemIdx = 0;
                    // fill result with computed gradient values : R
                    for (int i = 0; i < sz_.n_; ++i)
                    {
                        // ROFL_VAR1(gElemIdx);
                        // ROFL_VAR2("\n", rgR[gElemIdx]);
                        // result->GetElement(gElemIdx).SetToIdentity(); // Ri
                        // result->GetElement(gElemIdx).Print("Ri before assignment");

                        Vector RnextROPT(staircaseNextStepLevel, sz_.d_);
                        // RnextROPT.Initialize();
                        realdp *GroptlibWriteArray = RnextROPT.ObtainWriteEntireData();
                        for (int j = 0; j < rotSz; ++j)
                        {
                            // ROFL_VAR2(i, j);
                            // RnextROPT.Print("RnextROPT before assignment");

                            // ROFL_VAR1(RnextROPT.GetElement(j, 0));

                            GroptlibWriteArray[j] = Rnext[i].reshaped(sz_.d_ * staircaseNextStepLevel, 1)(j);

                            // ROFL_VAR1("");
                            // RnextROPT.Print("RnextROPT after assignment");
                        }
                        RnextROPT.CopyTo(xIn.GetElement(gElemIdx));
                        // result->GetElement(gElemIdx).Print("Riem. grad Ri after assignment");
                        gElemIdx++;
                    }

                    // fill result with computed gradient values : T

                    Vector TnextROPT(staircaseNextStepLevel, sz_.n_);
                    realdp *GroptlibWriteArray = TnextROPT.ObtainWriteEntireData();
                    for (int j = 0; j < staircaseNextStepLevel * sz_.n_; ++j)
                    {
                        // TnextROPT.Print("TnextROPT before assignment");

                        // ROFL_VAR1(RnextROPT.GetElement(j, 0));

                        GroptlibWriteArray[j] = Tnext.reshaped(sz_.n_ * staircaseNextStepLevel, 1)(j);

                        // ROFL_VAR1("");
                        // TnextROPT.Print("TnextROPT after assignment");
                    }
                    TnextROPT.CopyTo(xIn.GetElement(gElemIdx));
                } // end of EigToRopt scope for xIn

                linesearchArmijoROPTLIB(xIn, szNext, Y0);
            }
            else
            {
                SomUtils::MatD Y0T(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
                SomUtils::VecMatD Y0R(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
                linesearchDummy(Rnext, Tnext, vPimR, vPimT, Y0R, Y0T);
                Y0 = ProdManiNext.RandominManifold();
                // ROFL_VAR1(Y0T);

                { // EigToRopt scope for Y0

                    int rotSz = szNext.p_ * szNext.d_;
                    // int translSz = szNext.p_;

                    int gElemIdx = 0;
                    // fill result with computed gradient values : R
                    for (int i = 0; i < sz_.n_; ++i)
                    {
                        // ROFL_VAR1(gElemIdx);
                        // ROFL_VAR2("\n", rgR[gElemIdx]);
                        // result->GetElement(gElemIdx).SetToIdentity(); // Ri
                        // result->GetElement(gElemIdx).Print("Ri before assignment");

                        Vector RnextROPT(staircaseNextStepLevel, sz_.d_);
                        // RnextROPT.Initialize();
                        realdp *GroptlibWriteArray = RnextROPT.ObtainWriteEntireData();
                        for (int j = 0; j < rotSz; ++j)
                        {
                            // ROFL_VAR2(i, j);
                            // RnextROPT.Print("RnextROPT before assignment");

                            // ROFL_VAR1(RnextROPT.GetElement(j, 0));

                            GroptlibWriteArray[j] = Y0R[i].reshaped(sz_.d_ * staircaseNextStepLevel, 1)(j);

                            // ROFL_VAR1("");
                            // RnextROPT.Print("RnextROPT after assignment");
                        }
                        RnextROPT.CopyTo(Y0.GetElement(gElemIdx));
                        // result->GetElement(gElemIdx).Print("Riem. grad Ri after assignment");
                        gElemIdx++;
                    }

                    // fill result with computed gradient values : T

                    Vector TnextROPT(staircaseNextStepLevel, sz_.n_);
                    realdp *GroptlibWriteArray = TnextROPT.ObtainWriteEntireData();
                    for (int j = 0; j < staircaseNextStepLevel * sz_.n_; ++j)
                    {
                        // TnextROPT.Print("TnextROPT before assignment");

                        // ROFL_VAR1(RnextROPT.GetElement(j, 0));

                        GroptlibWriteArray[j] = Y0T.reshaped(sz_.n_ * staircaseNextStepLevel, 1)(j);

                        // ROFL_VAR1("");
                        // TnextROPT.Print("TnextROPT after assignment");
                    }
                    TnextROPT.Print("line 1215");
                    Y0.GetElement(gElemIdx).Print("line 1216");
                    TnextROPT.CopyTo(Y0.GetElement(gElemIdx));
                } // end of EigToRopt scope for xIn
            }
        }

        void linesearchArmijoROPTLIB(const Vector &xIn, const SomUtils::SomSize &somSzLocal, Vector &Y0) const
        {
            ROFL_VAR1("Running linesearchArmijoROPTLIB()")
            Stiefel mani1(somSzLocal.p_, somSzLocal.d_);
            mani1.ChooseParamsSet2();
            integer numoftypes = 2; // 2 i.e. (3D) Stiefel + Euclidean

            integer numofmani1 = somSzLocal.n_; // num of Stiefel manifolds
            integer numofmani2 = 1;
            Euclidean mani2(somSzLocal.p_, somSzLocal.n_);
            ProductManifold ProdMani(numoftypes, &mani1, numofmani1, &mani2, numofmani2);

            SampleSomProblem ProbLS(somSzLocal, Tijs_, edges_);
            ProbLS.SetDomain(&ProdMani);

            std::cout << "cost of LS input " << ProbLS.f(xIn) << std::endl; // x cost

            xIn.Print("Printing xIn inside linesearchArmijoROPTLIB()");

            ROPTLIB::RNewton *RNewtonSolver = new RNewton(&ProbLS, &xIn); // USE INITGUESS HERE!
            RNewtonSolver->Verbose = ITERRESULT;                          // TODO: comment this to remove some output

            PARAMSMAP solverParams = {std::pair<std::string, double>("Max_Inner_Iter", 1)};
            RNewtonSolver->SetParams(solverParams);
            RNewtonSolver->CheckParams();

            RNewtonSolver->LineSearch_LS = LSSM_ARMIJO;
            RNewtonSolver->LinesearchInput = &SomUtils::LinesearchInput;
            // LinesearchInput(iter, x1, eta1, stepsize, initialslope, Prob, this)
            // RNewtonSolver->IsPureLSInput = false;
            // RNewtonSolver->StopPtr = &MyStop;
            // RNewtonSolver->Max_Iteration = 25; // should be the same as the default one

            RNewtonSolver->Run();

            RNewtonSolver->GetXopt().Print("Xopt");

            Y0 = RNewtonSolver->GetXopt();

            std::cout << "cost of LS output " << RNewtonSolver->Getfinalfun() << std::endl; // x cost
        }

        void stiefelRetractionQR(const SomUtils::VecMatD &x, const SomUtils::VecMatD &e, SomUtils::VecMatD &rxe) const
        {
            // N = size(x, 3);
            // p = size(x, 2);
            // if N > 1
            //     rxe = zeros(size(x));
            //     for ii = 1 : N
            //         x_ii = x( :, :, ii);
            //         e_ii = e( :, :, ii);
            //         [ Q, R ] = qr(x_ii + e_ii);
            //         rxe( :, :, ii) = R;
            //     end
            // else % rxe = zeros(size(x));
            //     [ Q, R ] = qr(x + e);
            //     rxe = R;
            // end
        }

        void euclRetraction(const SomUtils::MatD &x, const SomUtils::MatD &d, SomUtils::MatD &y, double t = 1.0) const
        {
            ROFL_ASSERT(y.rows() == x.rows() && y.rows() == d.rows())
            ROFL_ASSERT(y.cols() == x.cols() && y.cols() == d.cols())
            y = x + t * d;
        }

        void stiefelRetraction(const SomUtils::MatD &xIn, const SomUtils::MatD &e, SomUtils::MatD &rxe) const
        {

            // else %     rxe = zeros(size(x));
            //     i_p = eye(p,p);
            //     snd_term = inv(sqrtm((i_p + e' * e)));
            //     rxe = (x + e)*snd_term;
            // end

            rxe.setZero();
            SomUtils::MatD Ip(SomUtils::MatD::Identity(sz_.p_, sz_.p_));
            SomUtils::MatD sndTerm = (Ip + e.transpose() * e).sqrt().inverse();
            rxe = (xIn + e) * sndTerm;

            // Assert check that new element is still on Stiefel
            ROFL_ASSERT(((rxe.transpose() * rxe) - SomUtils::MatD::Identity(sz_.d_, sz_.d_)).cwiseAbs().maxCoeff() < 1e-6);
        }

        void stiefelRetraction(const SomUtils::VecMatD &xIn, const SomUtils::VecMatD &e, SomUtils::VecMatD &rxe) const
        {
            // N = size(x, 3);
            // p = size(x, 2);

            // if N> 1
            //     rxe = zeros(size(x));
            //     i_p = eye3d(p, p, N);
            //     for ii = 1 : N
            //         x_ii = x( :, :, ii);
            //         e_ii = e( :, :, ii);
            //         snd_term_ii = inv(sqrtm((i_p(:,:,ii) + e_ii' * e_ii)));
            //         rxe(:,:,ii) = (x_ii + e_ii)*snd_term_ii;
            //     end

            for (int i = 0; i < sz_.n_; ++i)
            {
                stiefelRetraction(xIn[i], e[i], rxe[i]);
            }
        }

        void linesearchDummy(const SomUtils::VecMatD &xRin, const SomUtils::MatD &xTin,
                             const SomUtils::VecMatD &vRin, const SomUtils::MatD &vTin,
                             SomUtils::VecMatD &Y0R, SomUtils::MatD &Y0T) const
        {
            ROFL_VAR1("Running linesearchDummy()")

            int nrs = xTin.rows();
            SomUtils::SomSize szNext(nrs, sz_.d_, sz_.n_);

            stiefelRetraction(xRin, vRin, Y0R);
            ROFL_VAR1(Y0R[0]);

            euclRetraction(xTin, vTin, Y0T);
            ROFL_VAR1(Y0T)
        }

        // Code for recovery

        void makeAdjMatFromEdges(Eigen::MatrixXi &adjMat) const
        {
            ROFL_ASSERT(adjMat.rows() == numEdges_)
            ROFL_ASSERT(adjMat.cols() == numEdges_)

            adjMat.setZero();
            for (int k = 0; k < numEdges_; ++k)
            {
                int ii = edges_(k, 0) - 1;
                int jj = edges_(k, 1) - 1;
                adjMat(ii, jj) = 1;
            }
        }

        void computeNodeDegrees(Eigen::ArrayXi &nodeDegrees) const
        {
            ROFL_ASSERT(nodeDegrees.rows() == sz_.n_)

            Eigen::MatrixXi adjMat(Eigen::MatrixXi::Zero(numEdges_, numEdges_));
            makeAdjMatFromEdges(adjMat);

            nodeDegrees.setZero();
            nodeDegrees = adjMat.colwise().sum();
        }

        void makeTedges(const SomUtils::MatD &T, SomUtils::MatD &Tedges, SomUtils::MatD &T1offset) const
        {
            int nrs = T.rows();
            Tedges.setZero();
            ROFL_ASSERT(Tedges.rows() == nrs && Tedges.cols() == numEdges_)
            for (int e = 0; e < numEdges_; ++e)
            {
                int ii = edges_(e, 0) - 1;
                int jj = edges_(e, 1) - 1;

                Tedges.col(e) = T.col(ii) - T.col(jj);
            }

            ROFL_ASSERT(T1offset.rows() == nrs && T1offset.cols() == 1)
            T1offset = T.col(0);
        }

        void makeTedges(const SomUtils::MatD &T, SomUtils::MatD &Tedges) const
        {
            int nrs = T.rows();
            Tedges.setZero();
            ROFL_ASSERT(Tedges.rows() == nrs && Tedges.cols() == numEdges_)
            for (int e = 0; e < numEdges_; ++e)
            {
                int ii = edges_(e, 0) - 1;
                int jj = edges_(e, 1) - 1;

                Tedges.col(e) = T.col(ii) - T.col(jj);
            }
        }

        void POCRotateToMinimizeLastEntries(const SomUtils::MatD &x, SomUtils::MatD &Qtransp) const
        {
            Eigen::JacobiSVD<SomUtils::MatD> svd(x, Eigen::ComputeFullV | Eigen::ComputeFullU); // TODO: ComputeFullV flag can probably be removed
            auto Q = svd.matrixU();
            Qtransp = Q.transpose();
        }

        void makeTij1j2sEdges(int nodeId, const Eigen::ArrayXi &nodeDegrees, const SomUtils::MatD &Tedges,
                              SomUtils::MatD &Tij1j2, SomUtils::MatD &Tij1j2_tilde) const
        {
            // num_rows_T = size(T_edges, 1);
            int numRowsT = Tedges.rows();
            // % nrs = params.nrs;
            // d = params.d;
            auto nodeDeg = nodeDegrees(nodeId); // usually equal to low deg

            ROFL_ASSERT(Tij1j2.rows() == sz_.d_ && Tij1j2.cols() == nodeDeg)
            // SomUtils::MatD Tij1j2 = SomUtils::MatD::Zero(sz_.d_, nodeDeg);
            ROFL_ASSERT(Tij1j2_tilde.rows() == numRowsT && Tij1j2_tilde.cols() == nodeDeg)
            // SomUtils::MatD Tij1j2_tilde = SomUtils::MatD::Zero(num_rows_T, node_deg);

            Tij1j2.setZero();
            Tij1j2_tilde.setZero();

            int found = 1;
            for (int e = 0; e < edges_.rows(); ++e)
            {
                auto eI = edges_(e, 0);
                // % e_j = edges(e, 2);
                if (eI == nodeId)
                {
                    Tij1j2.col(found) = Tijs_.col(e);
                    Tij1j2_tilde.col(found) = -Tedges.col(e);
                    found++;
                }
            }
            ROFL_ASSERT(found == numEdges_)
        }

        void recoverRiTilde(const SomUtils::MatD &Ritilde2, const SomUtils::MatD &Tij1j2_tilde,
                            SomUtils::MatD &RitildeEst1, SomUtils::MatD &RitildeEst2) const
        {
            // PARAMS IN: Qx_edges*R_i_tilde2, Tij1j2_tilde

            // Qx = align2d(Tijtilde);
            // QxRitilde2Bot = Qx(3 : 4, :) * Ritilde2;
            // [ U, ~, ~] = svd(QxRitilde2Bot, 'econ');
            // c = U( :, 2);

            // QLastRight = Qx(3 : 4, 4)';

            //     RbEst = procrustesRb(c, QLastRight'); RitildeEst1 = Qx '*blkdiag(eye(2),-RbEst') * Qx * Ritilde2;
            // RitildeEst2=Qx'*blkdiag(eye(2),RbEst')*Qx*Ritilde2;
        }

        struct Vertex
        {
            int index;

            Vertex() : index(-1) {}

            Vertex(int idx) : index(idx) {}

            ~Vertex() {}
        };

        void dijkstraBT(int src, int n, const std::vector<int> &prev, std::vector<std::vector<int>> &list) const
        {
            list.clear();
            list.resize(n);

            for (int i = 0; i < n; ++i)
            {
                list[i].push_back(src);

                if (i == src)
                {
                    continue;
                }
                int curr = prev[i];
                int currIdx = i;
                do
                {
                    // ROFL_VAR3(i, curr, prev[curr]);
                    if (curr != src)
                        list[i].insert(list[i].begin() + 1, curr);
                    currIdx = curr;
                    curr = prev[curr];
                } while (currIdx != src);

                // list[i].pop_back();

                list[i].push_back(i);
            }
        }

        void dijkstraBTedges(int src, int n, const std::vector<int> &prev, const Eigen::MatrixXi &edges, std::vector<std::vector<int>> &listEdges) const
        {
            std::vector<std::vector<int>> listNodes;
            dijkstraBT(src, n, prev, listNodes);

            listEdges.clear();
            listEdges.resize(n);

            // for each node
            for (int i = 0; i < n; ++i)
            {
                if (i == src)
                    continue;

                auto listNodeI = listNodes[i];

                // run through node list (of i-th node)
                for (int j = 0; j < listNodeI.size() - 1; ++j)
                {
                    int idxI = listNodeI[j] + 1;     //+1!!
                    int idxJ = listNodeI[j + 1] + 1; //+1!!

                    // add edges
                    bool addedEdge = false;
                    for (int k = 0; k < edges.rows(); ++k)
                    {
                        if (edges(k, 0) == idxI && edges(k, 1) == idxJ)
                        {
                            listEdges[i].push_back(k);
                            addedEdge = true;
                            break;
                        }
                    }
                    ROFL_ASSERT(addedEdge)
                }
            }
        }

        void dijkstraSP(int n, int src, const Eigen::MatrixXi &adjmat, std::vector<double> &dist, std::vector<int> &prev) const
        {
            // function Dijkstra(Graph, source):

            // for each vertex v in Graph.Vertices:
            //     dist[v] ← INFINITY
            //     prev[v] ← UNDEFINED
            //     add v to Q

            dist.clear();
            prev.clear();
            dist.assign(n, INFINITY);
            prev.assign(n, nan("nan"));
            std::vector<Vertex> q(n);
            for (int i = 0; i < n; ++i)
            {
                Vertex qI(i);
                q[i] = qI;
            }

            // ROFL_VAR1("dist, prev init")
            // for (int i = 0; i < dist.size(); ++i)
            // {
            //     ROFL_VAR2(dist[i], prev[i]);
            // }

            // dist[source] ← 0

            dist[src] = 0.0;

            // while Q is not empty:
            while (!q.empty())
            {
                // ROFL_VAR1(q.size());

                //     u ← vertex in Q with minimum dist[u]
                Vertex uVert;
                double mindist = INFINITY;
                int u = -1;
                int iToErase = -1;
                for (int i = 0; i < q.size(); ++i)
                {
                    auto qI = q[i];
                    // ROFL_VAR2(i, qI.index)
                    if (dist[qI.index] < mindist)
                    {
                        iToErase = i;
                        // ROFL_VAR1("Saving i, qI.index")
                        uVert = qI;
                        u = qI.index;
                        mindist = dist[qI.index];
                    }
                }

                // ROFL_VAR2(u, mindist)

                //     remove u from Q
                q.erase(q.begin() + iToErase);

                // ROFL_VAR2("After erase", q.size())

                //     for each neighbor v of u still in Q:
                for (int v = 0; v < adjmat.cols(); ++v)
                {
                    bool isNeighborInQ = false;
                    for (int j = 0; j < q.size(); ++j)
                    {
                        if (adjmat(u, v) > 0 && q[j].index == v)
                        {
                            isNeighborInQ = true;
                            break;
                        }
                    }
                    if (!isNeighborInQ)
                        continue;

                    // ROFL_VAR2(q.size(), v);

                    //         alt ← dist[u] + Graph.Edges(u, v)
                    double alt = mindist + adjmat(u, v);
                    //         if alt < dist[v]:
                    if (alt < dist[v])
                    {
                        //             dist[v] ← alt
                        //             prev[v] ← u
                        // ROFL_VAR1("Updating dist, prev")
                        dist[v] = alt;
                        prev[v] = u;
                    }
                }

                // for (int i = 0; i < n; ++i)
                // {
                //     ROFL_VAR3(i, dist[i], prev[i]);
                // }
            }

            // return dist[], prev[]
        }

        void edgeDiffs2T(const SomUtils::MatD &Tdiffs, int n, SomUtils::MatD &T) const
        {
            // nrs = size(T_diffs, 1);
            int nrs = Tdiffs.rows();
            ROFL_ASSERT(nrs == T.rows())
            ROFL_ASSERT(T.cols() == n)
            // booleans_T = boolean(0) * ones(N,1); % alg should stop when all these are 1
            Eigen::ArrayXi booleansT(Eigen::ArrayXi::Ones(n));
            // booleans_T(1) = boolean(1); % node 1 chosen as reference
            booleansT(0) = true;
            // T = zeros(nrs, N);
            T.setZero();
            // adjmat = edges2adjmatrix(edges);
            Eigen::MatrixXi adjMat(Eigen::MatrixXi::Zero(numEdges_, numEdges_));
            makeAdjMatFromEdges(adjMat);
            //
            int globalId = 0;
            std::vector<double> dist;
            std::vector<int> prev;
            dijkstraSP(sz_.n_, globalId, adjMat, dist, prev);
            std::vector<std::vector<int>> listEdges;
            dijkstraBTedges(globalId, sz_.n_, prev, edges_, listEdges);

            for (int i = 0; i < sz_.n_; ++i)
            {
                if (i == globalId)
                    continue;
                // %     fprintf("shortest_p for node %g\n", ii);
                // %     disp(shortest_p)
                // %     fprintf("length for node %g\n", ii);
                // %     disp(length)
                // %     fprintf("edge_path for node %g\n", ii);
                // %     disp(edge_path)
                auto edgePath = listEdges[i];
                for (int j = 0; j < edgePath.size(); ++j)
                {
                    T.col(i) += Tdiffs.col(j);
                }
                T.col(i) *= -1;
                //     for ep = edge_path
                //         T(:,ii) = T(:,ii) + T_diffs(:,ep);
                //         disp('');
                //     end
                //     T(:,ii) = -T(:,ii);
            }
        }

        bool recoverySEdN(int staircaseStepIdx,
                          const SomUtils::VecMatD &RmanoptOut, const SomUtils::MatD &TmanoptOut,
                          SomUtils::VecMatD &Rrecovered, SomUtils::MatD &Trecovered) const
        {
            ROFL_ASSERT(Rrecovered.size() == sz_.n_ && Trecovered.rows() == sz_.d_ && Trecovered.cols() == sz_.n_)
            if (staircaseStepIdx > sz_.d_ + 1)
            {
                int nrs = staircaseStepIdx - 1;
                int lowDeg = 2; // TODO: not necessarily 2 in more complex graph cases (?)

                Eigen::ArrayXi nodeDegrees(Eigen::ArrayXi::Zero(sz_.n_));
                computeNodeDegrees(nodeDegrees);

                auto nodesHighDeg = nodeDegrees > lowDeg;

                SomUtils::MatD Tedges(SomUtils::MatD::Zero(nrs, numEdges_));

                // RT_stacked_high_deg = [ matStackH(R_manopt_out( :, :, nodes_high_deg)), T_edges ];
                auto numNodesHighDeg = nodesHighDeg.sum();
                SomUtils::MatD RTstackedHighDeg(SomUtils::MatD::Zero(nrs, sz_.d_ * numNodesHighDeg + Tedges.cols()));
                SomUtils::MatD RstackedHighDeg(SomUtils::MatD::Zero(nrs, sz_.d_ * numNodesHighDeg));
                hstack(RmanoptOut, RstackedHighDeg);
                RTstackedHighDeg.block(0, 0, nrs, sz_.d_ * numNodesHighDeg) = RstackedHighDeg;
                RTstackedHighDeg.block(0, sz_.d_ * numNodesHighDeg, nrs, numEdges_) = Tedges;

                SomUtils::MatD QxEdges(SomUtils::MatD::Zero(nrs, nrs));
                POCRotateToMinimizeLastEntries(RTstackedHighDeg, QxEdges);

                // R_tilde2_edges = multiprod(repmat(Qx_edges, 1, 1, sum(nodes_high_deg)), R_manopt_out( :, :, nodes_high_deg));
                SomUtils::VecMatD Rtilde2edges(numNodesHighDeg, SomUtils::MatD::Zero(nrs, sz_.d_));
                int highDegId = 0;
                for (int i = 0; i < sz_.n_; ++i)
                {
                    if (nodesHighDeg[i])
                    {
                        Rtilde2edges[highDegId] = QxEdges * RmanoptOut[i];
                        highDegId++;
                    }
                }
                ROFL_ASSERT(highDegId == numNodesHighDeg)

                std::for_each(Rrecovered.begin(), Rrecovered.end(), [](SomUtils::MatD &x) { //^^^ take argument by reference: LAMBDA FUNCTION
                    x.setZero();
                });

                // R_recovered( :, :, nodes_high_deg) = R_tilde2_edges(1 : d, :, :);
                highDegId = 0;
                for (int i = 0; i < sz_.n_; ++i)
                {
                    if (nodesHighDeg[i])
                    {
                        ROFL_ASSERT(Rrecovered[i].rows() == sz_.d_ && Rrecovered[i].cols() == sz_.d_)
                        Rrecovered[i] = Rtilde2edges[highDegId].block(0, 0, sz_.d_, sz_.d_);
                        highDegId++;
                    }
                }
                ROFL_ASSERT(highDegId == numNodesHighDeg)

                auto nodesLowDeg = !nodesHighDeg;
                ROFL_VAR1(nodesLowDeg.transpose())

                if (!nodesLowDeg.any())
                {
                    ROFL_VAR1("No nodes low deg!");
                    SomUtils::MatD TdiffsShifted = QxEdges * Tedges; // this has last row to 0
                    edgeDiffs2T(TdiffsShifted.block(0, 0, sz_.d_, TdiffsShifted.cols()), sz_.n_, Trecovered);
                }
                else
                {
                    // for node_id = 1 : length(params.node_degrees)
                    //
                    for (int nodeId = 0; nodeId < nodeDegrees.size(); ++nodeId)
                    {
                        // node_deg = params.node_degrees(node_id);
                        auto nodeDeg = nodeDegrees(nodeId);
                        // if node_deg == low_deg
                        if (nodeDeg == lowDeg)
                        {
                            //             fprintf("Running recoverRitilde() on node %g\n", node_id);
                            std::cout << "Running recoverRitilde() on node " << nodeId << std::endl;
                            // R_i_tilde2 = R_manopt_out( :, :, node_id);
                            auto RiTilde2 = RmanoptOut[nodeId];

                            SomUtils::MatD Xgt(SomUtils::MatD::Zero(sz_.d_, sz_.d_ * sz_.n_ + sz_.d_ * sz_.n_));
                            SomUtils::MatD RgtSt(SomUtils::MatD::Zero(sz_.d_, sz_.d_ * sz_.n_));
                            // ROFL_VAR1("hstack call from here");
                            hstack(Rgt_, RgtSt);
                            ROFL_VAR1(RgtSt);

                            Xgt.block(0, 0, sz_.d_, RgtSt.cols()) = RgtSt;
                            Xgt.block(0, RgtSt.cols(), sz_.d_, Tgt_.cols()) = Tgt_;

                            // disp("cost_gt")
                            // disp(cost_gt)
                            double costGt = costEigen(Xgt);
                            ROFL_VAR1(costGt)
                            

                            SomUtils::MatD XmanoptOut(SomUtils::MatD::Zero(nrs, sz_.d_ * sz_.n_ + sz_.d_ * sz_.n_));
                            SomUtils::MatD RmanoptOutSt(SomUtils::MatD::Zero(nrs, sz_.d_ * sz_.n_));
                            // ROFL_VAR1("hstack call from here");
                            hstack(RmanoptOut, RmanoptOutSt);
                            ROFL_VAR1(RmanoptOutSt);

                            XmanoptOut.block(0, 0, nrs, RmanoptOutSt.cols()) = RmanoptOutSt;
                            XmanoptOut.block(0, RmanoptOutSt.cols(), nrs, TmanoptOut.cols()) = TmanoptOut;

                            double costManoptOutput = costEigen(XmanoptOut); // TOCHECK: is costEigen actually defined for this kind of input?

                            // disp("cost_manopt_output")
                            // disp(cost_manopt_output)
                            ROFL_VAR1(costManoptOutput);
                            // T_diffs_shifted = Qx_edges * T_edges;
                            auto TdiffsShifted = QxEdges * Tedges; // this has last row to 0

                            // [~, Tij1j2_tilde] = make_Tij1j2s_edges(node_id, T_diffs_shifted, Tijs, edges, params);
                            SomUtils::MatD Tij1j2(SomUtils::MatD::Zero(sz_.d_, nodeDeg));
                            SomUtils::MatD Tij1j2_tilde(SomUtils::MatD::Zero(nrs, nodeDeg));
                            makeTij1j2sEdges(nodeId, nodeDegrees, Tedges, Tij1j2, Tij1j2_tilde);

                            // [ RitildeEst1, RitildeEst2, ~, ~] = recoverRitilde(Qx_edges * R_i_tilde2, Tij1j2_tilde);
                            SomUtils::MatD RiTildeEst1(SomUtils::MatD::Zero(nrs, sz_.d_));
                            SomUtils::MatD RiTildeEst2(SomUtils::MatD::Zero(nrs, sz_.d_));
                            recoverRiTilde(QxEdges * RiTilde2, Tij1j2_tilde, RiTildeEst1, RiTildeEst2); // TODO: add possibility of returning "local" Qx

                            // disp('')
                            std::cout << std::endl; // TODO : how to decide between RitildeEst1, RitildeEst2 ? ? det_RitildeEst1 = det(RitildeEst1(1 : d, :));
                            // det_RitildeEst2 = det(RitildeEst2(1 : d, :));
                            auto detRiTildeEst1 = RiTildeEst1.block(0, 0, sz_.d_, sz_.d_).determinant();
                            auto detRiTildeEst2 = RiTildeEst2.block(0, 0, sz_.d_, sz_.d_).determinant();

                            // use_positive_det = boolean(1);
                            bool usePositiveDet = true;

                            // if (sum(multidet(R_tilde2_edges(1 : d, :, :))) < 0)
                            //     use_positive_det = boolean(0);
                            double tmp = 0.0;
                            for (int i = 0; i < Rtilde2edges.size(); ++i)
                            {
                                tmp += Rtilde2edges[i].block(0, 0, sz_.d_, sz_.d_).determinant();
                            }
                            if (tmp < 0)
                                usePositiveDet = false;

                            if (detRiTildeEst1 > 1 - 1e-5 && detRiTildeEst1 < 1 + 1e-5)
                            {
                                ROFL_ASSERT(Rrecovered[nodeId].rows() == sz_.d_ && Rrecovered[nodeId].cols() == sz_.d_)

                                //      if use_positive_det
                                //          R_recovered( :, :, node_id) = RitildeEst1(1 : d, :);
                                //      else
                                //          R_recovered( :, :, node_id) = RitildeEst2(1 : d, :);
                                if (usePositiveDet)
                                    Rrecovered[nodeId] = RiTildeEst1.block(0, 0, sz_.d_, sz_.d_);
                                else
                                    Rrecovered[nodeId] = RiTildeEst2.block(0, 0, sz_.d_, sz_.d_);
                            }
                            else if (detRiTildeEst2 > 1 - 1e-5 && detRiTildeEst2 < 1 + 1e-5)
                            {
                                ROFL_ASSERT(Rrecovered[nodeId].rows() == sz_.d_ && Rrecovered[nodeId].cols() == sz_.d_)

                                //     if use_positive_det
                                //          R_recovered( :, :, node_id) = RitildeEst2(1 : d, :);
                                //     else
                                //          R_recovered( :, :, node_id) = RitildeEst1(1 : d, :);
                                if (usePositiveDet)
                                    Rrecovered[nodeId] = RiTildeEst2.block(0, 0, sz_.d_, sz_.d_);
                                else
                                    Rrecovered[nodeId] = RiTildeEst1.block(0, 0, sz_.d_, sz_.d_);
                            }
                            else
                            {
                                // if ~params.noisy_test
                                // {
                                ROFL_VAR1("ERROR in recovery: Ritilde DETERMINANTS ~= +-1\n")
                                // rs_recovery_success = boolean(0); // maybe add possibility to return this also
                                // save('data/zerodet_ws.mat')
                                ROFL_ASSERT(0) // TODO: remove this line for noisy cases (required condition is too strong)
                                // }
                            }
                            // T_recovered = edge_diffs_2_T(T_diffs_shifted(1 : d, :), edges, N);
                            edgeDiffs2T(TdiffsShifted, sz_.n_, Trecovered);
                            std::cout << std::endl;
                        }
                    }
                }
            }
            else
            {
                // recovery is not actually performed but using the same variable names for simplicity
                Rrecovered = RmanoptOut;
                Trecovered = TmanoptOut;
            }

            // checking that cost has not changed during "recovery" X_recovered.T = T_recovered;
            // X_recovered.R = R_recovered;
            // cost_out = rsom_cost_base(X_recovered, problem_struct_next);
            // disp("cost_out")
            // disp(cost_out)
            SomUtils::MatD Xrecovered(SomUtils::MatD::Zero(sz_.d_, sz_.d_ * sz_.n_ + sz_.d_ * sz_.n_));
            SomUtils::MatD RrecoveredSt(SomUtils::MatD::Zero(sz_.d_, sz_.d_ * sz_.n_));
            // ROFL_VAR1("hstack call from here");
            hstack(Rrecovered, RrecoveredSt);
            ROFL_VAR1(RrecoveredSt);
            Xrecovered.block(0, 0, sz_.d_, RrecoveredSt.cols()) = RrecoveredSt;
            Xrecovered.block(0, RrecoveredSt.cols(), sz_.d_, Trecovered.cols()) = Trecovered;

            ROFL_VAR1("Before end of recoverySEdN(9)")
            //  disp("[matStackH(X_gt.R); matStackH(R_recovered)]");
            //  disp([matStackH(X_gt.R); matStackH(R_recovered)]);
            for (int i=0; i<sz_.n_; ++i)
            {
                ROFL_VAR3(i, Rgt_[i], Rrecovered[i])
            }
        }

        bool isEqualFloats(const SomUtils::MatD &a, const SomUtils::MatD &b, double thr = 1e-5) const
        {
            ROFL_ASSERT(a.rows() == b.rows() && a.cols() == b.cols())

            double val = (a - b).cwiseAbs().maxCoeff();

            if (val > thr)
                return false;

            return true;
        }

        bool isEqualFloats(const SomUtils::VecMatD &a, const SomUtils::VecMatD &b, double thr = 1e-5) const
        {
            int n = a.size();
            ROFL_ASSERT(n == b.size())
            bool retval = true;
            for (int i = 0; i < n; ++i)
            {
                if (!isEqualFloats(a, b, thr))
                    return false;
            }
            return true;
        }

        void globalize(int src, const SomUtils::VecMatD &Rsedn, const SomUtils::MatD &Tsedn,
                       SomUtils::VecMatD &Rout, SomUtils::MatD &Tout)
        {
            // R_recovered -> Rsedn
            // T_recovered -> Tsedn

            // GLOBALIZATION! -> probably put it in another function?
            // R_global = R_recovered(:,:,1) * X_gt.R(:,:,1)'; %!!
            auto Rglobal = Rsedn[src] * Rgt_[src].transpose();

            // code for making all rotations global at once
            // R_recovered_global = multiprod(repmat(R_global', 1, 1, N), R_recovered);
            SomUtils::VecMatD RrecoveredGlobal(sz_.n_, SomUtils::MatD::Zero(sz_.d_, sz_.d_));
            for (int i = 0; i < sz_.n_; ++i)
            {
                RrecoveredGlobal[i] = Rglobal.transpose() * Rsedn[i];
            }
            // disp("[matStackH(X_gt.R); matStackH(R_recovered_global)]");
            // disp([matStackH(X_gt.R); matStackH(R_recovered_global)]);
            for (int i = 0; i < sz_.n_; ++i)
            {
                ROFL_VAR3(i, Rgt_[i], RrecoveredGlobal[i])
            }

            // T_global = R_global * T_recovered(:,1) - X_gt.T(:,1); %!!
            auto Tglobal = Rglobal * Tsedn.col(src) - Tgt_.col(src);
            // code for making all translation global at once
            // disp("[X_gt.T; T_recovered]");

            // T_recovered_global = R_global' * T_recovered - T_global;
            // disp([X_gt.T; T_recovered_global]);
            auto TrecoveredGlobal = Rglobal.transpose() * Tsedn - Tglobal;
            ROFL_VAR2(Tgt_, TrecoveredGlobal)

            // Checking recovery success
            // for ii = 1:N
            rsRecoverySuccess_ = true;
            for (int i = 0; i < sz_.n_; ++i)
            {
                //     R_gt_i = X_gt.R(:,:,ii);
                auto RgtI = Rgt_[i];
                //     R_recov_i_global = R_recovered_global(:,:,ii); %GLOBAL!
                auto RrecovIglobal = RrecoveredGlobal[i];
                //     fprintf("ii %g\n", ii);
                ROFL_VAR1(i)
                //     % rotations
                //     disp("R_gt_i, R_recov_i_global");
                //     disp([R_gt_i, R_recov_i_global]);
                //     disp("is_equal_floats(R_gt_i, R_recov_i_global)")
                //     disp(is_equal_floats(R_gt_i, R_recov_i_global))
                ROFL_VAR2(RgtI, RrecovIglobal)
                //     if (~is_equal_floats(R_gt_i, R_recov_i_global))
                // %         error("rot found NOT equal")
                //         fprintf("ERROR in recovery: R_GLOBAL\n");
                //         rs_recovery_success = boolean(0);
                if (!isEqualFloats(RgtI, RrecovIglobal))
                {
                    ROFL_VAR1("ERROR in recovery: R_GLOBAL")
                    rsRecoverySuccess_ = false;
                    // ROFL_ASSERT(0)
                }
                //     % translations
                //     T_gt_i = X_gt.T(:,ii);
                auto TgtI = Tgt_.col(i);
                //     T_recov_i_global = T_recovered_global(:,ii);
                auto TrecovIglobal = TrecoveredGlobal.col(i);
                //     disp("[X_gt.T, T_recovered]");
                //     disp([T_gt_i, T_recov_i_global]);
                ROFL_VAR2(TgtI, TrecovIglobal);

                //     disp("is_equal_floats(T_gt_i, T_recov_i_global)")
                //     disp(is_equal_floats(T_gt_i, T_recov_i_global))
                ROFL_VAR1(isEqualFloats(TgtI, TrecovIglobal))
                //     if (~is_equal_floats(T_gt_i, T_recov_i_global))
                // %         error("transl found NOT equal")
                //         fprintf("ERROR in recovery: T_GLOBAL\n");
                //         rs_recovery_success = boolean(0);
                if (!isEqualFloats(TgtI, TrecovIglobal))
                {
                    ROFL_VAR1("ERROR in recovery: T_GLOBAL")
                    rsRecoverySuccess_ = false;
                    // ROFL_ASSERT(0)
                }
            }

            // fprintf("rs_recovery_success: %g\n", rs_recovery_success);
            ROFL_VAR1(rsRecoverySuccess_)

            // X_recovered_global.R = R_recovered_global;
            // X_recovered_global.T = T_recovered_global;
            // cost_out_global = rsom_cost_base(X_recovered_global, problem_struct_next);
            // disp("cost_out_global")
            // disp(cost_out_global)
            SomUtils::MatD Xout(SomUtils::MatD::Zero(sz_.d_, sz_.d_ * sz_.n_ + sz_.d_ * sz_.n_));
            SomUtils::MatD RoutSt(SomUtils::MatD::Zero(sz_.d_, sz_.d_ * sz_.n_));
            // ROFL_VAR1("hstack call from here");
            hstack(Rout, RoutSt);
            ROFL_VAR1(RoutSt);
            Xout.block(0, 0, sz_.d_, RoutSt.cols()) = RoutSt;
            Xout.block(0, RoutSt.cols(), sz_.d_, Tout.cols()) = Tout;
            ROFL_VAR1(costEigen(Xout))

            // transf_out = RT2G(R_recovered_global, T_recovered_global); %rsom_genproc() function output

            // DETERMINANTS CHECK
            std::vector<double> multidetRrecovered, multidetRrecoveredGlobal;
            // disp('multidet(R_recovered)')
            // disp(multidet(R_recovered))
            multidet(Rsedn, multidetRrecovered);
            for (int i = 0; i < sz_.n_; ++i)
                ROFL_VAR2(i, multidetRrecovered[i]);

            // disp('multidet(R_recovered_global)')
            // disp(multidet(R_recovered_global))
            multidet(Rout, multidetRrecoveredGlobal);
            for (int i = 0; i < sz_.n_; ++i)
                ROFL_VAR2(i, multidetRrecoveredGlobal[i]);
        }

        void multidet(const SomUtils::VecMatD &a3d, std::vector<double> &dets) const
        {
            int n = a3d.size();
            dets.clear();
            dets.assign(0.0, n);
            ROFL_ASSERT(n == dets.size())

            for (int i = 0; i < n; ++i)
                dets[i] = a3d[i].determinant();
        }

        SomUtils::VecMatD Rgt_;

        SomUtils::MatD Tgt_;

        bool rsRecoverySuccess_;
    };

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