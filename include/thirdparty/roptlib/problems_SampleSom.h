// #include "manifolds_Rotations.h"
#include "manifolds_Stiefel.h"
#include "manifolds_Euclidean.h"
#include "manifolds_MultiManifolds.h"

#include <vector>

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
            mOut = mIn;
            double normF = mIn.norm(); // TODO: maybe use .normalized() directly?

            mOut /= normF;
            ROFL_VAR3(mIn, normF, mOut);
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
            tmp.Print("tmp inside stiefelRandTgNormVector()");
            SomUtils::MatD tmpEigVec(SomUtils::MatD::Zero(r * c, 1));
            RoptToEigStiefel(tmp, tmpEigVec); // xEigen is supposed to be a vectorized matrix
            ROFL_VAR1(tmpEigVec);
            SomUtils::MatD tmpEig = tmpEigVec.reshaped<Eigen::RowMajor>(r, c);
            ROFL_VAR1(tmpEig);
            SomUtils::MatD tmpProj(SomUtils::MatD::Zero(r, c));
            stiefelTangentProj(mIn, tmpEig, tmpProj);
            ROFL_VAR1(tmpProj);

            normalizeEucl(tmpProj, mOut);
            ROFL_VAR1(mOut);

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

            int fullIdx = 0;
            for (int i = 0; i < sz_.n_; ++i)
            {
                for (int j = 0; j < sz_.d_; ++j)
                {
                    for (int k = 0; k < sz_.p_; ++k)
                    {
                        RvecOut(fullIdx, 0) = R[i](k, j);
                        fullIdx++;
                        // ROFL_VAR4(i, j, k, fullIdx);
                    }
                }
            }
        }

        void pimFunctionGenproc(const SomUtils::VecMatD &xR, const SomUtils::MatD &xT, const SomUtils::VecMatD &uR, const SomUtils::MatD &uT,
                                double &lambdaMax, SomUtils::VecMatD &uFullR, SomUtils::MatD &uFullT, double thresh = 1e-5) const
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
                hessGenprocEigen(xR, uRunstackedTmp, xT, uTcopy, uRunstackedOutTmp, uTout);

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

            hessGenprocEigen(xR, uRunstacked, xT, uTcopy / normRTmax, fxRunstackedOut, uTout1);

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
            unStackH(uFullRhSt, uFullR, sz_.d_);
            uFullT = uFullHst.block(0, sz_.d_ * sz_.n_, staircaseLevel, sz_.n_);

            // lambda_max = lambda_max / sum(stiefel_metric([], full_xmax( :), full_xmax( :)));
            lambdaMax /= uFullHst.norm();
        }

        void pimFunctionGenprocShifted(const SomUtils::VecMatD &xR, const SomUtils::MatD &xT, const SomUtils::VecMatD &uR, const SomUtils::MatD &uT, double mu,
                                       double &lambdaMax, SomUtils::VecMatD &uFullR, SomUtils::MatD &uFullT, double thresh = 1e-5) const
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
            unStackH(uFullRhSt, uFullR, sz_.d_);
            uFullT = uFullHst.block(0, sz_.d_ * sz_.n_, staircaseLevel, sz_.n_);

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
                double mu = 100 * lambdaPim;

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
                eigencheckHessianGenprocShifted(highestNormEigenval, Rnext, vPimRshift, Tnext, vPimTshift, mu);
                ROFL_VAR3(highestNormEigenval, vPimRshift[0], vPimTshift);
            }
            else
            {
                highestNormEigenval = lambdaPim;
            }
        }

        void rsomPimHessianGenproc(double thresh, const SomUtils::VecMatD &R, const SomUtils::MatD &T, Vector &Y0) const
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
            int staircaseNextStepLevel = R.size() + 1;
            SomUtils::VecMatD Rnext(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            catZeroRow3dArray(R, Rnext);
            SomUtils::MatD Tnext(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            catZeroRow(T, Tnext);

            // stiefel_normalize_han = @(x) x./ (norm(x(:))); //Note: this is basically eucl_normalize_han

            // u_start.R = stiefel_randTangentNormVector(Rnext);
            SomUtils::VecMatD RnextTg(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            stiefelRandTgNormVector(Rnext, RnextTg);
            // u_start.R = stiefel_normalize(Rnext, u_start.R);
            SomUtils::VecMatD RnextTgNorm(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            normalizeEucl(RnextTg, RnextTgNorm);

            // u_start.T = rand(size(Tnext));
            auto TnextTg = SomUtils::MatD::Random(staircaseNextStepLevel, sz_.d_);
            // u_start.T = stiefel_normalize_han(u_start.T);
            SomUtils::MatD TnextTgNorm(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
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

                auto TnextTgShift = SomUtils::MatD::Random(staircaseNextStepLevel, sz_.d_);
                SomUtils::MatD TnextTgNormShift(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
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
                highestNormEigenval = lambdaPimShift + mu;
                eigencheckHessianGenprocShifted(highestNormEigenval, Rnext, vPimRshift, Tnext, vPimTshift, mu);
            } // SMALL VERSION UP TO HERE!!

            //////!!!!/////!!!!
            // // %Preparing linesearch
            // nrs_next = problem_struct_next.sz(1);
            // d = problem_struct_next.sz(2);
            // N = problem_struct_next.sz(3);
            int nrsNext = sz_.p_ + 1;
            SomUtils::SomSize szNext(nrsNext, sz_.d_, sz_.n_);

            Stiefel mani1(szNext.p_, szNext.d_);
            mani1.ChooseParamsSet2();
            integer numoftypes = 2; // 2 i.e. (3D) Stiefel + Euclidean

            integer numofmani1 = szNext.n_; // num of Stiefel manifolds
            integer numofmani2 = 1;
            Euclidean mani2(szNext.p_, szNext.n_);
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

            { // EigToRopt scope for xIn

                int rotSz = getRotSz();
                int translSz = getTranslSz();

                int gElemIdx = 0;
                // fill result with computed gradient values : R
                for (int i = 0; i < sz_.n_; ++i)
                {
                    // ROFL_VAR1(gElemIdx);
                    // ROFL_VAR2("\n", rgR[gElemIdx]);
                    // result->GetElement(gElemIdx).SetToIdentity(); // Ri
                    // result->GetElement(gElemIdx).Print("Ri before assignment");

                    Vector RnextROPT(sz_.p_, sz_.d_);
                    // RnextROPT.Initialize();
                    realdp *GroptlibWriteArray = RnextROPT.ObtainWriteEntireData();
                    for (int j = 0; j < rotSz; ++j)
                    {
                        // ROFL_VAR2(i, j);
                        // RnextROPT.Print("RnextROPT before assignment");

                        // ROFL_VAR1(RnextROPT.GetElement(j, 0));

                        GroptlibWriteArray[j] = Rnext[i].reshaped(sz_.d_ * sz_.p_, 1)(j);

                        // ROFL_VAR1("");
                        // RnextROPT.Print("RnextROPT after assignment");
                    }
                    RnextROPT.CopyTo(xIn.GetElement(gElemIdx));
                    // result->GetElement(gElemIdx).Print("Riem. grad Ri after assignment");
                    gElemIdx++;
                }

                // fill result with computed gradient values : T

                Vector TnextROPT(sz_.p_, sz_.n_);
                realdp *GroptlibWriteArray = TnextROPT.ObtainWriteEntireData();
                for (int j = 0; j < sz_.p_ * sz_.n_; ++j)
                {
                    // TnextROPT.Print("TnextROPT before assignment");

                    // ROFL_VAR1(RnextROPT.GetElement(j, 0));

                    GroptlibWriteArray[j] = Tnext.reshaped(sz_.n_ * sz_.p_, 1)(j);

                    // ROFL_VAR1("");
                    // TnextROPT.Print("TnextROPT after assignment");
                }
                TnextROPT.CopyTo(xIn.GetElement(gElemIdx));
            } // end of EigToRopt scope for xIn

            // Vector Y0;
            linesearchArmijoROPTLIB(xIn, szNext, Y0);
        }

        void linesearchArmijoROPTLIB(const Vector &xIn, const SomUtils::SomSize &somSzLocal, Vector &Y0) const
        {
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
            RNewtonSolver->Verbose = ITERRESULT;

            PARAMSMAP solverParams = {std::pair<std::string, double>("Max_Inner_Iter", 1)};
            RNewtonSolver->SetParams(solverParams);
            RNewtonSolver->CheckParams();

            RNewtonSolver->LineSearch_LS = LSSM_INPUTFUN;
            RNewtonSolver->LinesearchInput = &SomUtils::LinesearchInput;
            // LinesearchInput(iter, x1, eta1, stepsize, initialslope, Prob, this)
            RNewtonSolver->IsPureLSInput = false;
            // RNewtonSolver->StopPtr = &MyStop;
            RNewtonSolver->Max_Iteration = 1;

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

        void linesearchDummy(const SomUtils::VecMatD &xRin, const SomUtils::MatD &xTin, SomUtils::VecMatD &xRout, SomUtils::MatD &xTout)
        {
            // rsomPimHessianGenproc(double thresh, const SomUtils::VecMatD &R, const SomUtils::MatD &T, Vector &Y0);

            int staircaseNextStepLevel = xRin.size() + 1;
            // SomUtils::VecMatD Rnext(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            // catZeroRow3dArray(xRin, Rnext);
            // SomUtils::MatD Tnext(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            // catZeroRow(xTin, Tnext);

            // stiefel_normalize_han = @(x) x./ (norm(x(:))); //Note: this is basically eucl_normalize_han

            // u_start.R = stiefel_randTangentNormVector(Rnext);
            SomUtils::VecMatD RnextTg(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            stiefelRandTgNormVector(xRin, RnextTg);
            // u_start.R = stiefel_normalize(Rnext, u_start.R);
            SomUtils::VecMatD RnextTgNorm(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            normalizeEucl(RnextTg, RnextTgNorm);

            // u_start.T = rand(size(Tnext));
            auto TnextTg = SomUtils::MatD::Random(staircaseNextStepLevel, sz_.d_);
            // u_start.T = stiefel_normalize_han(u_start.T);
            SomUtils::MatD TnextTgNorm(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            normalizeEucl(TnextTg, TnextTgNorm);

            // [lambda_pim, v_pim] = pim_function_genproc(rhess_fun_han, u_start, stiefel_normalize_han, thresh);
            // disp('Difference between lambda*v_max and H(v_max) should be in the order of the tolerance:')
            double lambdaPim = 1e+6;
            SomUtils::VecMatD vPimR(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            SomUtils::MatD vPimT(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));

            pimFunctionGenproc(xRin, xTin, RnextTgNorm, TnextTgNorm, lambdaPim, vPimR, vPimT);
            std::cout << "Difference between lambda_pim_after_shift*v_pim_after_shift"
                      << " and H_SH(v_pim_after_shift) should be in the order of the tolerance:" << std::endl;
            eigencheckHessianGenproc(lambdaPim, xRin, vPimR, xTin, vPimT);

            // if lambda_pim>0
            double highestNormEigenval = 1e+6;
            SomUtils::VecMatD vPimRshift(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            SomUtils::MatD vPimTshift(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
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
                stiefelRandTgNormVector(xRin, RnextTgShift);
                SomUtils::VecMatD RnextTgNormShift(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
                normalizeEucl(RnextTgShift, RnextTgNormShift);

                auto TnextTgShift = SomUtils::MatD::Random(staircaseNextStepLevel, sz_.d_);
                SomUtils::MatD TnextTgNormShift(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
                normalizeEucl(TnextTgShift, TnextTgNormShift);

                double lambdaPimShift = 1e+6; // "after" shift is intended
                pimFunctionGenprocShifted(xRin, xTin, RnextTgNorm, TnextTgNorm, mu, lambdaPimShift, vPimRshift, vPimTshift);

                //     disp(['Difference between lambda_pim_after_shift*v_pim_after_shift ' ...
                //         'and H_SH(v_pim_after_shift) should be in the order of the tolerance:'])
                //     eigencheck_hessian_genproc(lambda_pim_after_shift, v_pim_after_shift, ...
                //         rhess_shifted_fun_han);
                std::cout << "Difference between lambda_pim_after_shift*v_pim_after_shift"
                          << "and H_SH(v_pim_after_shift) should be in the order of the tolerance:" << std::endl;
                highestNormEigenval = lambdaPimShift + mu;
                eigencheckHessianGenprocShifted(highestNormEigenval, xRin, vPimRshift, xTin, vPimTshift, mu);
            }

            //////!!!!/////!!!!
            // // %Preparing linesearch
            // nrs_next = problem_struct_next.sz(1);
            // d = problem_struct_next.sz(2);
            // N = problem_struct_next.sz(3);
            int nrsNext = sz_.p_ + 1;
            SomUtils::SomSize szNext(nrsNext, sz_.d_, sz_.n_);

            SomUtils::VecMatD Y0R(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            stiefelRetraction(xRin, vPimRshift, Y0R);

            // Need also eucl. retraction
        };
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