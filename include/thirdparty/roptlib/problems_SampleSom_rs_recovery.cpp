#include "problems_SampleSom.h"

namespace ROPTLIB
{
    ////////////////////////////////////////RS////////////////////////////////////////
    bool SampleSomProblem::checkIsStiefelTg(const SomUtils::MatD &m, const SomUtils::MatD &mTg) const
    {
        // For any matrix representative $U \in St(n, p)$, the tangent space of $St(n, p)$ at $U$ is represented by
        // U\transpose \Delta = -\Delta\transpose U
        double diff = (m.transpose() * mTg + mTg.transpose() * m).cwiseAbs().maxCoeff();
        if (!(diff >= 0 && diff < 1e-5))
        {
            ROFL_VAR1("Check failed here")
            return false;
        }

        return true;
    }

    bool SampleSomProblem::checkIsStiefelTg(const SomUtils::VecMatD &m, const SomUtils::VecMatD &mTg) const
    {
        int n = m.size();
        ROFL_ASSERT(n == mTg.size())

        for (int i = 0; i < n; ++i)
        {
            if (!checkIsStiefelTg(m[i], mTg[i]))
                return false;
        }

        return true;
    }

    bool SampleSomProblem::checkIsStiefelTgNorm(const SomUtils::MatD &m, const SomUtils::MatD &mTg) const
    {
        // For any matrix representative $U \in St(n, p)$, the tangent space of $St(n, p)$ at $U$ is represented by
        // U\transpose \Delta = -\Delta\transpose U
        double diff = (m.transpose() * mTg + mTg.transpose() * m).cwiseAbs().maxCoeff();
        if (!(diff >= 0 && diff < 1e-5))
        {
            ROFL_VAR1("Check failed here")
            return false;
        }

        double mOutNorm = mTg.norm();
        if (!(mOutNorm > 1 - 1e-5 && mOutNorm < 1 + 1e-5, mOutNorm))
        {
            ROFL_VAR1("Check failed here")
            return false;
        }

        return true;
    }

    void SampleSomProblem::stiefelRandTgNormVector(const SomUtils::MatD &mIn, SomUtils::MatD &mOut) const
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
        SomUtils::stiefelTangentProj(mIn, tmpEig, tmpProj);
        // ROFL_VAR1(tmpProj);

        SomUtils::normalizeEucl(tmpProj, mOut);
        // ROFL_VAR1(mOut);

        ROFL_ASSERT(checkIsStiefelTg(mIn, mOut))
    }

    void SampleSomProblem::stiefelRandTgNormVector(const SomUtils::VecMatD &mIn, SomUtils::VecMatD &mOut) const
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

    bool SampleSomProblem::eigencheckHessianGenproc(const double &lambda,
                                                    const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR,
                                                    const SomUtils::MatD &xT, const SomUtils::MatD &uT, double thr) const
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

    bool SampleSomProblem::eigencheckHessianGenprocShifted(const double &lambda,
                                                           const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR,
                                                           const SomUtils::MatD &xT, const SomUtils::MatD &uT, double mu, double thr) const
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

    void SampleSomProblem::pimFunctionGenproc(const SomUtils::VecMatD &xR, const SomUtils::MatD &xT, const SomUtils::VecMatD &uR, const SomUtils::MatD &uT,
                                              double &lambdaMax, SomUtils::VecMatD &uOutR, SomUtils::MatD &uOutT, double thresh) const
    {
        ROFL_VAR1("Running pimFunctionGenproc()")

        // Note: normalization is done across entire ProdMani vector through simple eucl. metric

        // % % R iterative_change = 1e+6;
        // xR = x_start.R;
        // xT = x_start.T;
        // xfull = [ matStackH(x_start.R), x_start.T ];

        int staircaseLevel = xT.rows();
        ROFL_VAR1(staircaseLevel);

        SomUtils::MatD uRhStacked(SomUtils::MatD::Zero(staircaseLevel, sz_.d_ * sz_.n_));
        // ROFL_VAR1("hstack call from here");
        SomUtils::hstack(uR, uRhStacked);
        // ROFL_VAR1(uRhStacked);

        SomUtils::MatD uTcopy = uT; // useful for keeping const in function params
        // ROFL_VAR1(uTcopy);

        SomUtils::MatD uFullHst(SomUtils::MatD::Zero(staircaseLevel, sz_.n_ + uRhStacked.cols()));
        uFullHst.block(0, 0, staircaseLevel, uRhStacked.cols()) = uRhStacked;
        uFullHst.block(0, uRhStacked.cols(), staircaseLevel, sz_.n_) = uTcopy;
        // ROFL_VAR1(uFullHst);

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
        while (iterationNum < 5000) // && iterativeChange < 1e-3
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
            SomUtils::unStackH(uRhStacked, uRunstackedTmp);

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
            SomUtils::hstack(uRunstackedOutTmp, uRhStacked);
            // ROFL_VAR1("hstack call from here");
            uFullHst.block(0, 0, staircaseLevel, uRhStacked.cols()) = uRhStacked;
            uFullHst.block(0, uRhStacked.cols(), staircaseLevel, sz_.n_) = uTcopy;
            // ROFL_VAR1(uFullHst);

            //      iterative_change = max(normalization_fun(xfull_prev - xfull), [], "all");
            double iterativeChange = (uFullHstPrev - uFullHst).cwiseAbs().maxCoeff();
            // ROFL_VAR2(iterationNum, iterativeChange);
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
        SomUtils::unStackH(uRhStacked / normRTmax, uRunstacked, sz_.d_);
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
        SomUtils::unStackH(uFullRhSt, uOutR, sz_.d_);
        uOutT = uFullHst.block(0, sz_.d_ * sz_.n_, staircaseLevel, sz_.n_);

        // 6
        // lambda_max = lambda_max / sum(stiefel_metric([], full_xmax( :), full_xmax( :)));
        lambdaMax /= uFullHst.norm();
    }

    void SampleSomProblem::pimFunctionGenprocShifted(const SomUtils::VecMatD &xR, const SomUtils::MatD &xT, const SomUtils::VecMatD &uR, const SomUtils::MatD &uT, double mu,
                                                     double &lambdaMax, SomUtils::VecMatD &uOutR, SomUtils::MatD &uOutT, double thresh) const
    {
        ROFL_VAR1("Running pimFunctionGenprocShifted()")
        // Note: normalization is done across entire ProdMani vector through simple eucl. metric

        // % % R iterative_change = 1e+6;
        // xR = x_start.R;
        // xT = x_start.T;
        // xfull = [ matStackH(x_start.R), x_start.T ];

        int staircaseLevel = xT.rows();

        SomUtils::MatD uRhStacked(SomUtils::MatD::Zero(staircaseLevel, sz_.d_ * sz_.n_));
        // ROFL_VAR1("hstack call from here");
        SomUtils::hstack(uR, uRhStacked);

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
        while (iterationNum < 5000) // && iterativeChange < 1e-3
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
            SomUtils::unStackH(uRhStacked, uRunstackedTmp);

            SomUtils::VecMatD uRunstackedOutTmp(sz_.n_, SomUtils::MatD::Zero(staircaseLevel, sz_.d_));
            SomUtils::MatD uTout(SomUtils::MatD::Zero(staircaseLevel, sz_.n_));
            hessGenprocEigenShifted(xR, uRunstackedTmp, xT, uTcopy, mu, uRunstackedOutTmp, uTout);

            std::for_each(uRunstackedOutTmp.begin(), uRunstackedOutTmp.end(), [](SomUtils::MatD &x) { //^^^ take argument by reference: LAMBDA FUNCTION
                x *= -1;
            });
            uTcopy = -uTout;

            // ROFL_VAR1("hstack call from here");
            SomUtils::hstack(uRunstackedOutTmp, uRhStacked);
            uFullHst.block(0, 0, staircaseLevel, uRhStacked.cols()) = uRhStacked;
            uFullHst.block(0, uRhStacked.cols(), staircaseLevel, sz_.n_) = uTcopy;

            //      iterative_change = max(normalization_fun(xfull_prev - xfull), [], "all");
            double iterativeChange = (uFullHstPrev - uFullHst).cwiseAbs().maxCoeff();
            // ROFL_VAR2(iterationNum, iterativeChange);
        }

        // norm_RT_max = norm([ matStackH(x.R), x.T ]);
        double normRTmax = uFullHst.norm();

        // x_max.R = x.R / norm_RT_max;
        // x_max.T = x.T / norm_RT_max;
        uFullHst /= normRTmax;

        // f_x_max = f(x_max);
        SomUtils::VecMatD uRunstacked(sz_.n_, SomUtils::MatD::Zero(staircaseLevel, sz_.d_));
        SomUtils::unStackH(uRhStacked / normRTmax, uRunstacked, sz_.d_);
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
        SomUtils::unStackH(uFullRhSt, uOutR, sz_.d_);
        uOutT = uFullHst.block(0, sz_.d_ * sz_.n_, staircaseLevel, sz_.n_);

        // lambda_max = lambda_max / sum(stiefel_metric([], full_xmax( :), full_xmax( :)));
        lambdaMax /= uFullHst.norm();
    }

    void SampleSomProblem::rsomPimHessianGenprocSmall(double thresh, const SomUtils::VecMatD &R, const SomUtils::MatD &T) const
    {
        // // [Y_star, lambda, v] = rsom_pim_hessian_genproc( ...
        // //     X, problem_struct_next, thr);
        // // disp("v") // %just to remove unused variable warning
        // // disp(v)
        // // if lambda > 0
        // //     disp("R, T eigenvals > 0: exiting staircase")
        // // break;

        // /////////////////////////////////////////////////////
        // // if ~exist('thresh', 'var')
        // //     thresh = 1e-6;
        // // end

        // // Rnext = cat_zero_rows_3d_array(X.R);
        // // Tnext = cat_zero_row(X.T);
        // // Xnext.R = Rnext;
        // // Xnext.T = Tnext;
        // // rhess_fun_han = @(u) hess_genproc(Xnext,u,problem_struct_next);
        // int staircaseNextStepLevel = T.rows() + 1;
        // ROFL_VAR1(staircaseNextStepLevel);
        // SomUtils::VecMatD Rnext(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
        // catZeroRow3dArray(R, Rnext);
        // SomUtils::MatD Tnext(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
        // catZeroRow(T, Tnext);

        // // stiefel_normalize_han = @(x) x./ (norm(x(:))); //Note: this is basically eucl_normalize_han

        // // u_start.R = stiefel_randTangentNormVector(Rnext);
        // SomUtils::VecMatD RnextTg(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
        // stiefelRandTgNormVector(Rnext, RnextTg);
        // // u_start.R = stiefel_normalize(Rnext, u_start.R);
        // SomUtils::VecMatD RnextTgNorm(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
        // SomUtils::normalizeEucl(RnextTg, RnextTgNorm);

        // // u_start.T = rand(size(Tnext));
        // auto TnextTg = SomUtils::MatD::Random(staircaseNextStepLevel, sz_.n_);
        // // u_start.T = stiefel_normalize_han(u_start.T);
        // SomUtils::MatD TnextTgNorm(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
        // SomUtils::normalizeEucl(TnextTg, TnextTgNorm);

        // // [lambda_pim, v_pim] = pim_function_genproc(rhess_fun_han, u_start, stiefel_normalize_han, thresh);
        // // disp('Difference between lambda*v_max and H(v_max) should be in the order of the tolerance:')
        // double lambdaPim = 1e+6;
        // SomUtils::VecMatD vPimR(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
        // SomUtils::MatD vPimT(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));

        // pimFunctionGenproc(Rnext, Tnext, RnextTgNorm, TnextTgNorm, lambdaPim, vPimR, vPimT);
        // std::cout << "Difference between lambda_pim_after_shift*v_pim_after_shift"
        //           << " and H_SH(v_pim_after_shift) should be in the order of the tolerance:" << std::endl;
        // eigencheckHessianGenproc(lambdaPim, Rnext, vPimR, Tnext, vPimT);

        // // if lambda_pim>0
        // double highestNormEigenval = 1e+6;
        // if (lambdaPim > 0)
        // {
        //     std::cout << "lambdaPim " << lambdaPim << std::endl;
        //     double mu = 1.1 * lambdaPim;

        //     //     rhess_shifted_fun_han = @(u) hess_genproc_shifted(Xnext,u,mu,problem_struct_next);

        //     //     // %run shifted power iteration
        //     //     u_start_second_iter.R = stiefel_randTangentNormVector(Rnext);
        //     //     u_start_second_iter.R = stiefel_normalize(Rnext, u_start_second_iter.R);
        //     //     u_start_second_iter.T = rand(size(Tnext));
        //     //     u_start_second_iter.T = stiefel_normalize_han(u_start.T);
        //     //     [lambda_pim_after_shift, v_pim_after_shift] = pim_function_genproc( ...
        //     //         rhess_shifted_fun_han, u_start_second_iter, stiefel_normalize_han, thresh);
        //     SomUtils::VecMatD RnextTgShift(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
        //     stiefelRandTgNormVector(Rnext, RnextTgShift);
        //     SomUtils::VecMatD RnextTgNormShift(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
        //     SomUtils::normalizeEucl(RnextTgShift, RnextTgNormShift);

        //     auto TnextTgShift = SomUtils::MatD::Random(staircaseNextStepLevel, sz_.n_);
        //     SomUtils::MatD TnextTgNormShift(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
        //     SomUtils::normalizeEucl(TnextTgShift, TnextTgNormShift);

        //     SomUtils::VecMatD vPimRshift(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
        //     SomUtils::MatD vPimTshift(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
        //     double lambdaPimShift = 1e+6; // "after" shift is intended
        //     pimFunctionGenprocShifted(Rnext, Tnext, RnextTgNorm, TnextTgNorm, mu, lambdaPimShift, vPimRshift, vPimTshift);

        //     //     disp(['Difference between lambda_pim_after_shift*v_pim_after_shift ' ...
        //     //         'and H_SH(v_pim_after_shift) should be in the order of the tolerance:'])
        //     //     eigencheck_hessian_genproc(lambda_pim_after_shift, v_pim_after_shift, ...
        //     //         rhess_shifted_fun_han);
        //     std::cout << "Difference between lambda_pim_after_shift*v_pim_after_shift"
        //               << "and H_SH(v_pim_after_shift) should be in the order of the tolerance:" << std::endl;
        //     eigencheckHessianGenprocShifted(lambdaPimShift, Rnext, vPimRshift, Tnext, vPimTshift, mu);
        //     highestNormEigenval = lambdaPimShift + mu;
        //     std::cout << "Difference between (lambda_pim_after_shift + mu)*v_pim_after_shift"
        //               << "and H_SH(v_pim_after_shift) should be in the order of the tolerance:" << std::endl;
        //     eigencheckHessianGenproc(highestNormEigenval, Rnext, vPimRshift, Tnext, vPimTshift, mu);
        //     // ROFL_VAR3(highestNormEigenval, vPimRshift[0], vPimTshift);
        //     vPimR = vPimRshift;
        //     vPimT = vPimTshift;
        // }
        // else
        // {
        //     highestNormEigenval = lambdaPim;
        // }
    }

    void SampleSomProblem::rsomPimHessianGenprocEigen(double thresh,
                                                      const SomUtils::VecMatD &R, const SomUtils::MatD &T,
                                                      Vector &Y0, double &lambdaPimOut, SomUtils::VecMatD &vPimRout, SomUtils::MatD &vPimTout,
                                                      bool armijo) const
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
        SomUtils::catZeroRow3dArray(R, Rnext);
        SomUtils::MatD Tnext(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
        SomUtils::catZeroRow(T, Tnext);

        // stiefel_normalize_han = @(x) x./ (norm(x(:))); //Note: this is basically eucl_normalize_han

        // u_start.R = stiefel_randTangentNormVector(Rnext);
        SomUtils::VecMatD RnextTg(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
        stiefelRandTgNormVector(Rnext, RnextTg);
        // u_start.R = stiefel_normalize(Rnext, u_start.R);
        SomUtils::VecMatD RnextTgNorm(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
        SomUtils::normalizeEucl(RnextTg, RnextTgNorm);

        // u_start.T = rand(size(Tnext));
        auto TnextTg = SomUtils::MatD::Random(staircaseNextStepLevel, sz_.n_);
        // u_start.T = stiefel_normalize_han(u_start.T);
        SomUtils::MatD TnextTgNorm(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
        SomUtils::normalizeEucl(TnextTg, TnextTgNorm);

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
            stiefelRandTgNormVector(Rnext, RnextTgShift);
            SomUtils::VecMatD RnextTgNormShift(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
            SomUtils::normalizeEucl(RnextTgShift, RnextTgNormShift);

            auto TnextTgShift = SomUtils::MatD::Random(staircaseNextStepLevel, sz_.n_);
            SomUtils::MatD TnextTgNormShift(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
            SomUtils::normalizeEucl(TnextTgShift, TnextTgNormShift);

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
            // vPimR = vPimRshift;
            // vPimT = vPimTshift;
        }
        else
        {
            highestNormEigenval = lambdaPim;
            vPimRshift = vPimR;
            vPimTshift = vPimT;
        } // SMALL VERSION UP TO HERE!!

        if (highestNormEigenval > 0)
        {
            lambdaPimOut = highestNormEigenval;
            vPimRout = vPimRshift;
            vPimTout = vPimTshift;
            ROFL_VAR2(highestNormEigenval, "hne > 0 -> avoiding linesearch")
            return;
        }

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
            linesearchDummy(costCurr_, Rnext, Tnext, vPimRshift, vPimTshift, Y0R, Y0T);

            // ls Dummy debug: save on file scope
            {
                { // costCurr
                    std::ofstream outfile;

                    outfile.open("costCurr.csv", std::ios_base::trunc);
                    outfile << costCurr_;
                    outfile.close();
                }

                { // Rnext
                    std::ofstream outfile;

                    outfile.open("Rnext.csv", std::ios_base::trunc);
                    if (!outfile.is_open())
                    {
                        ROFL_VAR1("Error opening outfile")
                        ROFL_ASSERT(0);
                    }
                    for (int i = 0; i < sz_.n_; ++i)
                        outfile << Rnext[i].reshaped<Eigen::ColMajor>(Rnext[i].rows() * Rnext[i].cols(), 1) << std::endl;
                    outfile.close();
                }

                { // Tnext
                    std::ofstream outfile;

                    outfile.open("Tnext.csv", std::ios_base::trunc);
                    if (!outfile.is_open())
                    {
                        ROFL_VAR1("Error opening outfile")
                        ROFL_ASSERT(0);
                    }
                    outfile << Tnext.reshaped<Eigen::ColMajor>(Tnext.rows() * Tnext.cols(), 1);
                    outfile.close();
                }

                { // vPimRshift
                    std::ofstream outfile;

                    outfile.open("vPimRshift.csv", std::ios_base::trunc);
                    if (!outfile.is_open())
                    {
                        ROFL_VAR1("Error opening outfile")
                        ROFL_ASSERT(0);
                    }
                    for (int i = 0; i < sz_.n_; ++i)
                        outfile << vPimRshift[i].reshaped<Eigen::ColMajor>(vPimRshift[i].rows() * vPimRshift[i].cols(), 1) << std::endl;
                    outfile.close();
                }

                { // vPimTshift
                    std::ofstream outfile;

                    outfile.open("vPimTshift.csv", std::ios_base::trunc);
                    if (!outfile.is_open())
                    {
                        ROFL_VAR1("Error opening outfile")
                        ROFL_ASSERT(0);
                    }
                    outfile << vPimTshift.reshaped<Eigen::ColMajor>(vPimTshift.rows() * vPimTshift.cols(), 1);
                    outfile.close();
                }

                { // Y0R
                    std::ofstream outfile;

                    outfile.open("Y0R.csv", std::ios_base::trunc);
                    if (!outfile.is_open())
                    {
                        ROFL_VAR1("Error opening outfile")
                        ROFL_ASSERT(0);
                    }
                    for (int i = 0; i < sz_.n_; ++i)
                        outfile << Y0R[i].reshaped<Eigen::ColMajor>(Y0R[i].rows() * Y0R[i].cols(), 1) << std::endl;
                    outfile.close();
                }

                { // Y0T
                    std::ofstream outfile;

                    outfile.open("Y0T.csv", std::ios_base::trunc);
                    if (!outfile.is_open())
                    {
                        ROFL_VAR1("Error opening outfile")
                        ROFL_ASSERT(0);
                    }
                    outfile << Y0T.reshaped<Eigen::ColMajor>(Y0T.rows() * Y0T.cols(), 1);
                    outfile.close();
                }
            }

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
                // TnextROPT.Print("line 1215");
                // Y0.GetElement(gElemIdx).Print("line 1216");
                TnextROPT.CopyTo(Y0.GetElement(gElemIdx));
            } // end of EigToRopt scope for xIn
        }

        lambdaPimOut = highestNormEigenval;
        vPimRout = vPimRshift;
        vPimTout = vPimTshift;
    }

    void SampleSomProblem::rsomEscapeHessianGenprocEigen(const SomUtils::VecMatD &R, const SomUtils::MatD &T,
                                                         Vector &Y0, double &lambdaOut, SomUtils::VecMatD &vRout, SomUtils::MatD &vTout) const
    {
        int staircaseStepLevel = T.rows() + 1; // TODO: this has likely to be passed as param
        ROFL_VAR1(staircaseStepLevel);
        SomUtils::VecMatD Rnext(sz_.n_, SomUtils::MatD::Zero(staircaseStepLevel, sz_.d_));
        SomUtils::catZeroRow3dArray(R, Rnext);
        SomUtils::MatD Tnext(SomUtils::MatD::Zero(staircaseStepLevel, sz_.n_));
        SomUtils::catZeroRow(T, Tnext);

        /**
         * makeHmat()
         */

        // Xvec = vectorizeXrt(X);
        // vecsz = length(Xvec);
        // Hmat = zeros(vecsz);
        SomUtils::MatD Xvec(staircaseStepLevel * sz_.d_ * sz_.n_ + staircaseStepLevel * sz_.n_, 1);
        vectorizeRT(Rnext, Tnext, Xvec);

        int vecsz = Xvec.rows();
        SomUtils::MatD Hmat(SomUtils::MatD::Zero(vecsz, vecsz));

        SomUtils::SomSize somSzNext(staircaseStepLevel, sz_.d_, sz_.n_);

        ROPTLIB::SampleSomProblem ProbNext(somSzNext, Tijs_, edges_);

        ProbNext.makeHmat(Xvec, somSzNext, Hmat);
        ROFL_VAR1(Hmat)

        if (!SomUtils::isEqualFloats(Hmat - Hmat.transpose(), SomUtils::MatD::Zero(vecsz, vecsz)))
        {
            // Checking whether anti-symmetric part is zero
            ROFL_ERR("Hessian NOT symmetric")
            // ROFL_ASSERT(0)
        }

        // Eigen::SelfAdjointEigenSolver<SomUtils::MatD> es;
        // es.compute(Hmat);
        // auto eigvals = es.eigenvalues(); // TODO: check how to use sparse matrix methods -> maybe SPECTRA library?
        // auto eigvecs = es.eigenvectors();
        // ROFL_VAR5(vecsz, eigvals.rows(), eigvals.cols(), eigvecs.rows(), eigvecs.cols());
        // ROFL_VAR1(Hmat)

        /** Computing min eigencouple via Eigen library*/
        // ROFL_VAR1("Computing Eigenvalues and Eigenvectors");
        // Eigen::EigenSolver<SomUtils::MatD> es;
        // es.compute(Hmat);
        // auto eigvals = es.eigenvalues(); // TODO: check how to use sparse matrix methods -> maybe SPECTRA library?
        // auto eigvecs = es.eigenvectors();

        // ROFL_VAR2(eigvals, eigvecs);

        // Eigen::Index minRow, minCol;
        // lambdaOut = eigvals.real().minCoeff(&minRow, &minCol);

        // ROFL_VAR1(lambdaOut);

        // auto vOut = eigvecs.real().col(minRow);
        // ROFL_VAR2(vOut.rows(), vOut.cols());
        /** */

        /** Computing min eigencouple via Spectra library*/
        Eigen::SparseMatrix<double, Eigen::ColMajor> HmatSparse(Hmat.rows(), Hmat.cols());
        HmatSparse = Hmat.sparseView();
        // Eigen::EigenSolver<Eigen::SparseMatrix<double>> es;
        // es.compute(HmatSparse);
        // auto eigvals = es.eigenvalues();
        // auto eigvecs = es.eigenvectors();

        // Construct matrix operation object using the wrapper class DenseSymMatProd
        Spectra::SparseGenMatProd<double> op(HmatSparse);

        // Construct eigen solver object, requesting the largest three eigenvalues
        Spectra::GenEigsSolver<Spectra::SparseGenMatProd<double>> eigs(op, 1, HmatSparse.rows());

        // Initialize and compute
        eigs.init();
        int nconv = eigs.compute(Spectra::SortRule::SmallestReal);

        // Retrieve results
        Eigen::VectorXd eigvals;
        Eigen::MatrixXd eigvecs;
        if (eigs.info() == Spectra::CompInfo::Successful)
        {
            eigvals = eigs.eigenvalues().real();
            eigvecs = eigs.eigenvectors().real();
        }

        std::cout << "Eigenvalues found:\n"
                  << eigvals << std::endl;
        std::cout << "Eigenvectors found:\n"
                  << eigvecs << std::endl;

        auto vOut = eigvecs.real().col(0);        

        /** */

        lambdaOut = eigvals(0);
        ProbNext.getRotations(vOut, vRout);
        ProbNext.getTranslations(vOut, vTout);

        ProbNext.eigencheckHessianGenproc(lambdaOut, Rnext, vRout, Tnext, vTout);

        // auto minusVrOut = vRout;
        // std::for_each(minusVrOut.begin(), minusVrOut.end(), [](SomUtils::MatD &x) { //^^^ take argument by reference: LAMBDA FUNCTION
        //     x *= -1;
        // });
        // ProbNext.eigencheckHessianGenproc(lambdaOut, Rnext, minusVrOut, Tnext, -vTout);

        if (lambdaOut > -1e-2)
        {
            // lambdaOut, vRout, vTout are already set before this check
            ROFL_VAR2(lambdaOut, "lambdaOut > 0 -> avoiding linesearch")
            return;
        }

        ROFL_VAR4(eigvals.rows(), eigvals.cols(), eigvecs.rows(), eigvecs.cols());

        /**
         * Linesearch for escaping saddle point
         */

        // // %Preparing linesearch
        // nrs_next = problem_struct_next.sz(1);
        // d = problem_struct_next.sz(2);
        // N = problem_struct_next.sz(3);

        Stiefel mani1(staircaseStepLevel, somSzNext.d_);
        mani1.ChooseParamsSet2();
        integer numoftypes = 2; // 2 i.e. (3D) Stiefel + Euclidean

        integer numofmani1 = somSzNext.n_; // num of Stiefel manifolds
        integer numofmani2 = 1;
        Euclidean mani2(staircaseStepLevel, somSzNext.n_);
        ProductManifold ProdManiNext(numoftypes, &mani1, numofmani1, &mani2, numofmani2);

        Vector xIn = ProdManiNext.RandominManifold(); //!! in other cases xIn would have been a pointer

        SomUtils::MatD Y0T(SomUtils::MatD::Zero(staircaseStepLevel, sz_.n_));
        SomUtils::VecMatD Y0R(sz_.n_, SomUtils::MatD::Zero(staircaseStepLevel, sz_.d_));
        linesearchDummy(costCurr_, Rnext, Tnext, vRout, vTout, Y0R, Y0T);

        // ls Dummy debug: save on file scope
        {
            { // costCurr
                std::ofstream outfile;

                outfile.open("costCurr.csv", std::ios_base::trunc);
                outfile << costCurr_;
                outfile.close();
            }

            { // Rnext
                std::ofstream outfile;

                outfile.open("Rnext.csv", std::ios_base::trunc);
                if (!outfile.is_open())
                {
                    ROFL_VAR1("Error opening outfile")
                    ROFL_ASSERT(0);
                }
                for (int i = 0; i < sz_.n_; ++i)
                    outfile << Rnext[i].reshaped<Eigen::ColMajor>(Rnext[i].rows() * Rnext[i].cols(), 1) << std::endl;
                outfile.close();
            }

            { // Tnext
                std::ofstream outfile;

                outfile.open("Tnext.csv", std::ios_base::trunc);
                if (!outfile.is_open())
                {
                    ROFL_VAR1("Error opening outfile")
                    ROFL_ASSERT(0);
                }
                outfile << Tnext.reshaped<Eigen::ColMajor>(Tnext.rows() * Tnext.cols(), 1);
                outfile.close();
            }

            { // vRout
                std::ofstream outfile;

                outfile.open("vRout.csv", std::ios_base::trunc);
                if (!outfile.is_open())
                {
                    ROFL_VAR1("Error opening outfile")
                    ROFL_ASSERT(0);
                }
                for (int i = 0; i < sz_.n_; ++i)
                    outfile << vRout[i].reshaped<Eigen::ColMajor>(vRout[i].rows() * vRout[i].cols(), 1) << std::endl;
                outfile.close();
            }

            { // vTout
                std::ofstream outfile;

                outfile.open("vTout.csv", std::ios_base::trunc);
                if (!outfile.is_open())
                {
                    ROFL_VAR1("Error opening outfile")
                    ROFL_ASSERT(0);
                }
                outfile << vTout.reshaped<Eigen::ColMajor>(vTout.rows() * vTout.cols(), 1);
                outfile.close();
            }

            { // Y0R
                std::ofstream outfile;

                outfile.open("Y0R.csv", std::ios_base::trunc);
                if (!outfile.is_open())
                {
                    ROFL_VAR1("Error opening outfile")
                    ROFL_ASSERT(0);
                }
                for (int i = 0; i < sz_.n_; ++i)
                    outfile << Y0R[i].reshaped<Eigen::ColMajor>(Y0R[i].rows() * Y0R[i].cols(), 1) << std::endl;
                outfile.close();
            }

            { // Y0T
                std::ofstream outfile;

                outfile.open("Y0T.csv", std::ios_base::trunc);
                if (!outfile.is_open())
                {
                    ROFL_VAR1("Error opening outfile")
                    ROFL_ASSERT(0);
                }
                outfile << Y0T.reshaped<Eigen::ColMajor>(Y0T.rows() * Y0T.cols(), 1);
                outfile.close();
            }
        }

        Y0 = ProdManiNext.RandominManifold();
        // ROFL_VAR1(Y0T);

        { // EigToRopt scope for Y0

            int rotSz = somSzNext.p_ * somSzNext.d_;
            // int translSz = somSzNext.p_;

            int gElemIdx = 0;
            // fill result with computed gradient values : R
            for (int i = 0; i < sz_.n_; ++i)
            {
                // ROFL_VAR1(gElemIdx);
                // ROFL_VAR2("\n", rgR[gElemIdx]);
                // result->GetElement(gElemIdx).SetToIdentity(); // Ri
                // result->GetElement(gElemIdx).Print("Ri before assignment");

                Vector RnextROPT(staircaseStepLevel, sz_.d_);
                // RnextROPT.Initialize();
                realdp *GroptlibWriteArray = RnextROPT.ObtainWriteEntireData();
                for (int j = 0; j < rotSz; ++j)
                {
                    // ROFL_VAR2(i, j);
                    // RnextROPT.Print("RnextROPT before assignment");

                    // ROFL_VAR1(RnextROPT.GetElement(j, 0));

                    GroptlibWriteArray[j] = Y0R[i].reshaped(sz_.d_ * staircaseStepLevel, 1)(j);

                    // ROFL_VAR1("");
                    // RnextROPT.Print("RnextROPT after assignment");
                }
                RnextROPT.CopyTo(Y0.GetElement(gElemIdx));
                // result->GetElement(gElemIdx).Print("Riem. grad Ri after assignment");
                gElemIdx++;
            }

            // fill result with computed gradient values : T

            Vector TnextROPT(staircaseStepLevel, sz_.n_);
            realdp *GroptlibWriteArray = TnextROPT.ObtainWriteEntireData();
            for (int j = 0; j < staircaseStepLevel * sz_.n_; ++j)
            {
                // TnextROPT.Print("TnextROPT before assignment");

                // ROFL_VAR1(RnextROPT.GetElement(j, 0));

                GroptlibWriteArray[j] = Y0T.reshaped(sz_.n_ * staircaseStepLevel, 1)(j);

                // ROFL_VAR1("");
                // TnextROPT.Print("TnextROPT after assignment");
            }
            // TnextROPT.Print("line 1215");
            // Y0.GetElement(gElemIdx).Print("line 1216");
            TnextROPT.CopyTo(Y0.GetElement(gElemIdx));
        } // end of EigToRopt scope for xIn
    }

    void SampleSomProblem::rsomPimHessianGenproc(double thresh, const SomUtils::VecMatD &R, const SomUtils::MatD &T, Vector &Y0, bool armijo) const
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
        SomUtils::catZeroRow3dArray(R, Rnext);
        SomUtils::MatD Tnext(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
        SomUtils::catZeroRow(T, Tnext);

        // stiefel_normalize_han = @(x) x./ (norm(x(:))); //Note: this is basically eucl_normalize_han

        // u_start.R = stiefel_randTangentNormVector(Rnext);
        SomUtils::VecMatD RnextTg(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
        stiefelRandTgNormVector(Rnext, RnextTg);
        // u_start.R = stiefel_normalize(Rnext, u_start.R);
        SomUtils::VecMatD RnextTgNorm(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
        SomUtils::normalizeEucl(RnextTg, RnextTgNorm);

        // u_start.T = rand(size(Tnext));
        auto TnextTg = SomUtils::MatD::Random(staircaseNextStepLevel, sz_.n_);
        // u_start.T = stiefel_normalize_han(u_start.T);
        SomUtils::MatD TnextTgNorm(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
        SomUtils::normalizeEucl(TnextTg, TnextTgNorm);

        // [lambda_pim, v_pim] = pim_function_genproc(rhess_fun_han, u_start, stiefel_normalize_han, thresh);
        // disp('Difference between lambda*v_max and H(v_max) should be in the order of the tolerance:')
        double lambdaPim = 1e+6;
        SomUtils::VecMatD vPimR(sz_.n_, SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.d_));
        SomUtils::MatD vPimT(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));

        pimFunctionGenproc(Rnext, Tnext, RnextTgNorm, TnextTgNorm, lambdaPim, vPimR, vPimT);
        std::cout << "Difference between lambda_pim_after_shift*v_pim_after_shift"
                  << " and H_SH(v_pim_after_shift) should be in the order of the tolerance:" << std::endl;
        eigencheckHessianGenproc(lambdaPim, Rnext, vPimR, Tnext, vPimT);
        ROFL_VAR1(lambdaPim)

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
            SomUtils::normalizeEucl(RnextTgShift, RnextTgNormShift);

            auto TnextTgShift = SomUtils::MatD::Random(staircaseNextStepLevel, sz_.n_);
            SomUtils::MatD TnextTgNormShift(SomUtils::MatD::Zero(staircaseNextStepLevel, sz_.n_));
            SomUtils::normalizeEucl(TnextTgShift, TnextTgNormShift);

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
            ROFL_VAR1(lambdaPimShift)
            highestNormEigenval = lambdaPimShift + mu;
            std::cout << "Difference between (lambda_pim_after_shift + mu)*v_pim_after_shift"
                      << "and H_SH(v_pim_after_shift) should be in the order of the tolerance:" << std::endl;
            eigencheckHessianGenproc(highestNormEigenval, Rnext, vPimRshift, Tnext, vPimTshift, mu);
            ROFL_VAR1(highestNormEigenval)
            // ROFL_VAR3(highestNormEigenval, vPimRshift[0], vPimTshift);
            vPimR = vPimRshift;
            vPimT = vPimTshift;
            for (int i = 0; i < sz_.n_; ++i)
                ROFL_VAR2(i, vPimR[i])
            ROFL_VAR1(vPimT)
        }
        else
        {
            highestNormEigenval = lambdaPim;
        } // SMALL VERSION UP TO HERE!!

        if (highestNormEigenval > 0)
        {
            // lambdaPimOut = highestNormEigenval;
            // vPimRout = vPimRshift;
            // vPimTout = vPimTshift;
            ROFL_VAR2(highestNormEigenval, "hne > 0 -> avoiding linesearch")
            return;
        }

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
            linesearchDummy(costCurr_, Rnext, Tnext, vPimR, vPimT, Y0R, Y0T);
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
                // TnextROPT.Print("line 1215");
                // Y0.GetElement(gElemIdx).Print("line 1216");
                TnextROPT.CopyTo(Y0.GetElement(gElemIdx));
            } // end of EigToRopt scope for xIn
        }
    }

    void SampleSomProblem::linesearchArmijoROPTLIB(const Vector &xIn, const SomUtils::SomSize &somSzLocal, Vector &Y0) const
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

    // void SampleSomProblem::stiefelRetractionQR(const SomUtils::VecMatD &x, const SomUtils::VecMatD &e, SomUtils::VecMatD &rxe) const
    // {
    //     // N = size(x, 3);
    //     // p = size(x, 2);
    //     // if N > 1
    //     //     rxe = zeros(size(x));
    //     //     for ii = 1 : N
    //     //         x_ii = x( :, :, ii);
    //     //         e_ii = e( :, :, ii);
    //     //         [ Q, R ] = qr(x_ii + e_ii);
    //     //         rxe( :, :, ii) = R;
    //     //     end
    //     // else % rxe = zeros(size(x));
    //     //     [ Q, R ] = qr(x + e);
    //     //     rxe = R;
    //     // end
    // }

    // void SampleSomProblem::QRunique(const SomUtils::MatD &a, SomUtils::MatD &Q, SomUtils::MatD &R) const
    // {
    //     // [q, r] = qr(A(:, :, k), 0);

    //     // % [q,r] = QR(a) performs a QR decomposition on m-by-n matrix a such that a = q*r.
    //     // The factor r is an m-by-n upper triangular matrix and
    //     // q is an m-by-m unitary matrix.

    //     R.setIdentity(); // TODO: add R computation

    //     int m = a.rows();
    //     int n = a.cols();

    //     Eigen::ColPivHouseholderQR<SomUtils::MatD> qr(a);
    //     SomUtils::MatD q = qr.matrixQ();
    //     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> r = qr.matrixQR().triangularView<Eigen::Upper>();

    //     // ROFL_ASSERT_VAR4(SomUtils::isEqualFloats(a, q * r), SomUtils::isEqualFloats(a, q * r), a, q, r);
    //     // ROFL_ASSERT(r.isUpperTriangular() && r.rows() == m && r.cols() == n);
    //     ROFL_ASSERT(q.isUnitary() && q.rows() == m && q.cols() == m);

    //     // % In the real case, s holds the signs of the diagonal entries of R.
    //     // % In the complex case, s holds the unit-modulus phases of these
    //     // % entries. In both cases, d = diag(s) is a unitary matrix, and
    //     // % its inverse is d* = diag(conj(s)).
    //     // s = sign(diag(r));

    //     // % Since a = qr (with 'a' the slice of A currently being processed),
    //     // % it is also true that a = (qd)(d*r). By construction, qd still has
    //     // % orthonormal columns, and d*r has positive real entries on its
    //     // % diagonal, /unless/ s contains zeros. The latter can only occur if
    //     // % slice a does not have full column rank, so that the decomposition
    //     // % is not unique: we make an arbitrary choice in that scenario.
    //     // % While exact zeros are unlikely, they may occur if, for example,
    //     // % the slice a contains repeated columns, or columns that are equal
    //     // % to zero. If an entry should be mathematically zero but is only
    //     // % close to zero numerically, then it is attributed an arbitrary
    //     // % sign dictated by the numerical noise: this is also fine.
    //     // s(s == 0) = 1;

    //     SomUtils::MatD rDiag = r.diagonal();
    //     SomUtils::MatD s = rDiag;
    //     for (int i = 0; i < rDiag.size(); ++i)
    //     {
    //         if (rDiag(i, 0) >= 0)
    //             s(i, 0) = 1;
    //         else
    //             s(i, 0) = -1;
    //     }

    //     // Q(:, :, k) = bsxfun(@times, q, s.');
    //     for (int i = 0; i < Q.rows(); ++i)
    //     {
    //         for (int j = 0; j < Q.cols(); ++j)
    //         {
    //             auto sTransp = s.transpose();
    //             Q(i, j) = q(i, j) * sTransp(0, j);
    //         }
    //     }

    //     // R(:, :, k) = bsxfun(@times, r, conj(s));
    //     for (int i = 0; i < R.rows(); ++i)
    //     {
    //         for (int j = 0; j < R.cols(); ++j)
    //         {
    //             // TODO: add complex s case
    //             // auto sConj = s.conjugate();
    //             R(i, j) = R(i, j) * s(i, 0);
    //         }
    //     }
    //     // TODO: add complex s case
    //     // ROFL_ASSERT(SomUtils::isEqualFloats(R.imag(), SomUtils::MatD::Zero(R.imag().rows(), R.imag().cols())))
    //     // R = R.real();
    // }

    void SampleSomProblem::QRunique(const SomUtils::VecMatD &A, SomUtils::VecMatD &Q, SomUtils::VecMatD &R) const
    {
        // [m, n, N] = size(A);
        // if m >= n % A (or its slices) has more rows than columns
        //     Q = zeros(m, n, N, class(A));
        //     R = zeros(n, n, N, class(A));
        // else
        //     Q = zeros(m, m, N, class(A));
        //     R = zeros(m, n, N, class(A));
        // end

        int m = A[0].rows();
        int n = A[0].cols();
        int N = A.size();
        if (m >= n)
        {
            Q.resize(N, SomUtils::MatD::Zero(m, n));
            R.resize(N, SomUtils::MatD::Zero(n, n));
        }
        else
        {
            Q.resize(N, SomUtils::MatD::Zero(m, m));
            R.resize(N, SomUtils::MatD::Zero(m, n));
        }

        // for k = 1 : N
        // ...

        for (int i = 0; i < N; ++i)
        {
            ROFL_ASSERT(A[i].rows() == m && A[i].cols() == n)
            main_qr_unique(A[i], Q[i], R[i]);
        }
    }

    bool SampleSomProblem::checkIsOn3dStiefel(const SomUtils::VecMatD &m) const
    {
        int p = m[0].rows();
        int d = m[0].cols();

        int n = m.size();

        for (int i = 0; i < n; ++i)
        {
            ROFL_ASSERT(m[i].rows() == p && m[i].cols() == d)

            // ROFL_VAR1("Calling SomUtils::isEqualFloats()")
            if (!checkIsOnStiefel(m[i]))
                return false;
        }
        return true;
    }

    bool SampleSomProblem::checkIsOnStiefel(const SomUtils::MatD &m) const
    {
        int p = m.rows();
        int d = m.cols();

        // ROFL_VAR1("Calling SomUtils::isEqualFloats()")

        return SomUtils::isEqualFloats(m.transpose() * m, SomUtils::MatD::Identity(d, d));
    }

    void SampleSomProblem::stiefelRetractionQR(const SomUtils::VecMatD &x, const SomUtils::VecMatD &e, SomUtils::VecMatD &rxe, double t) const
    {
        // function Y = retraction_qr(X, U, t)
        //     % It is necessary to call qr_unique rather than simply qr to ensure
        //     % this is a retraction, to avoid spurious column sign flips.
        //     if nargin < 3
        //         Y = qr_unique(X + U);
        //     else
        //         Y = qr_unique(X + t*U);
        //     end
        // end

        ROFL_ASSERT_VAR2(checkIsStiefelTg(x, e), x[0], e[0])

        SomUtils::VecMatD rTmp(sz_.n_);
        SomUtils::VecMatD xeSum(sz_.n_, SomUtils::MatD::Zero(x[0].rows(), x[0].cols()));
        SomUtils::VecMatD te = e;
        std::for_each(te.begin(), te.end(), [t](SomUtils::MatD &x) { //^^^ take argument by reference: LAMBDA FUNCTION
            // ROFL_VAR1(alpha)
            x *= t;
        });
        std::transform(x.begin(), x.end(), te.begin(), xeSum.begin(), std::plus<SomUtils::MatD>());
        QRunique(xeSum, rxe, rTmp); // rTmp is unused; rxe comes out from Q in QR decomposition and is the function's output
        //!! line above would need main_qr_unique() with sz_.n_ = 1

        // for (int i = 0; i < sz_.n_; ++i)
        // {
        //     ROFL_VAR3(i, xeSum[i], rxe[i])
        // }
        ROFL_ASSERT(checkIsOn3dStiefel(rxe))

        // for (int i = 0; i < sz_.n_; ++i)
        // {
        // ROFL_VAR2(i, "Equal to matlab stiefel retraction qr?")
        // std::cout << std::endl
        //           << "x[i]" << std::endl;
        // std::cout << std::endl
        //           << x[i] << std::endl;
        // std::cout << std::endl
        //           << "e[i]" << std::endl;
        // std::cout << std::endl
        //           << e[i] << std::endl;
        // std::cout << std::endl
        //           << "rxe[i]" << std::endl;
        // std::cout << std::endl
        //           << rxe[i] << std::endl;
        // }
    }

    void SampleSomProblem::euclRetraction(const SomUtils::MatD &x, const SomUtils::MatD &d, SomUtils::MatD &y, double t) const
    {
        ROFL_ASSERT(y.rows() == x.rows() && y.rows() == d.rows())
        ROFL_ASSERT(y.cols() == x.cols() && y.cols() == d.cols())
        y = x + t * d;
    }

    void SampleSomProblem::stiefelRetractionPolar(const SomUtils::MatD &xIn, const SomUtils::MatD &e, SomUtils::MatD &rxe, double t) const
    {

        // else %     rxe = zeros(size(x));
        //     i_p = eye(p,p);
        //     snd_term = inv(sqrtm((i_p + e' * e)));
        //     rxe = (x + e)*snd_term;
        // end

        ROFL_ASSERT(xIn.rows() == e.rows() && xIn.rows() == rxe.rows())
        ROFL_ASSERT(xIn.cols() == e.cols() && xIn.cols() == rxe.cols())

        // rxe.setZero();
        // SomUtils::MatD Ip(SomUtils::MatD::Identity(sz_.p_, sz_.p_));
        // ROFL_VAR3(Ip, e, e.transpose() * e);
        // SomUtils::MatD sndTerm = (Ip + e.transpose() * e).sqrt().inverse();
        // rxe = (xIn + e) * sndTerm;

        // [u, s, v] = svd(Y(:, :, kk), 'econ'); %#ok
        // Y(:, :, kk) = u*v';

        Eigen::JacobiSVD<SomUtils::MatD> svd(xIn + e, Eigen::ComputeFullV | Eigen::ComputeThinU); // TODO: ComputeFullV flag can probably be removed
        auto U = svd.matrixU();
        auto V = svd.matrixV();
        // ROFL_VAR3(xIn + e, U, V)
        rxe = U * V.transpose();
        // Assert check that new element is still on Stiefel
        ROFL_ASSERT(((rxe.transpose() * rxe) - SomUtils::MatD::Identity(sz_.d_, sz_.d_)).cwiseAbs().maxCoeff() < 1e-6);
    }

    void SampleSomProblem::stiefelRetractionPolar(const SomUtils::VecMatD &xIn, const SomUtils::VecMatD &e, SomUtils::VecMatD &rxe, double t) const
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
            stiefelRetractionPolar(xIn[i], e[i], rxe[i]);
        }
    }

    void SampleSomProblem::linesearchDummy(const double costInit,
                                           const SomUtils::VecMatD &xRin, const SomUtils::MatD &xTin,
                                           const SomUtils::VecMatD &vRin, const SomUtils::MatD &vTin,
                                           SomUtils::VecMatD &Y0R, SomUtils::MatD &Y0T,
                                           bool qr) const
    {
        ROFL_VAR1("Running linesearchDummy()")

        int nrs = xTin.rows();
        SomUtils::SomSize szNext(nrs, sz_.d_, sz_.n_);

        Y0R = xRin;
        Y0T = xTin;

        SomUtils::VecMatD Y0Rtry = Y0R;
        SomUtils::MatD Y0Ttry = Y0T;

        int maxLsSteps = 25;            // TODO: add this as param
        double contractionFactor = 0.5; // TODO: add this as param

        SomUtils::MatD vVec(SomUtils::MatD::Zero(vRin[0].rows() * vRin[0].cols() * vRin.size() + vTin.rows() * vTin.cols(), 1));
        vectorizeRT(vRin, vTin, vVec);

        double dnorm = vVec.norm();

        ROFL_VAR1(dnorm);

        double initialStepsize = dnorm;
        double alpha = initialStepsize / dnorm;
        double f0 = costInit;

        SomUtils::VecMatD alphaVrIn = vRin;
        std::for_each(alphaVrIn.begin(), alphaVrIn.end(), [alpha](SomUtils::MatD &x) { //^^^ take argument by reference: LAMBDA FUNCTION
            // ROFL_VAR1(alpha)
            x *= alpha;
        });
        if (qr)
            stiefelRetractionQR(xRin, alphaVrIn, Y0Rtry);
        else
            stiefelRetractionPolar(xRin, alphaVrIn, Y0Rtry);
        // ROFL_VAR1(Y0R[0]);

        SomUtils::MatD alphaVtIn = alpha * vTin;

        euclRetraction(xTin, alphaVtIn, Y0Ttry);
        // ROFL_VAR1(Y0T)

        double newf = costEigen(Y0Rtry, Y0Ttry);
        ROFL_VAR2(newf, f0)

        int costEvaluations = 1;

        while (newf >= f0)
        {
            // Reduce the step size,
            alpha *= contractionFactor;

            // and look closer down the line.
            alphaVrIn = vRin;
            // ROFL_VAR1(alpha)
            std::for_each(alphaVrIn.begin(), alphaVrIn.end(), [alpha](SomUtils::MatD &x) { //^^^ take argument by reference: LAMBDA FUNCTION
                // ROFL_VAR1(alpha)
                x *= alpha;
            });
            alphaVtIn = alpha * vTin;

            if (qr)
                stiefelRetractionQR(xRin, alphaVrIn, Y0Rtry);
            else
                stiefelRetractionPolar(xRin, alphaVrIn, Y0Rtry);

            euclRetraction(xTin, alphaVtIn, Y0Ttry);

            newf = costEigen(Y0Rtry, Y0Ttry);

            ROFL_VAR2(newf, f0)

            costEvaluations++;

            // make sure we don't run out of budget.
            if (costEvaluations >= maxLsSteps)
                break;
        }

        if (newf <= f0)
        {
            // accept step
            Y0R = Y0Rtry;
            Y0T = Y0Ttry;
        }
    }

    ////////////////////////////////////////RECOVERY////////////////////////////////////////

    void SampleSomProblem::computeNodeDegrees(Eigen::ArrayXi &nodeDegrees) const
    {
        ROFL_ASSERT(nodeDegrees.rows() == sz_.n_)

        Eigen::MatrixXi adjMat(Eigen::MatrixXi::Zero(sz_.n_, sz_.n_));
        makeAdjMatFromEdges(adjMat);

        ROFL_VAR1(adjMat)

        nodeDegrees.setZero();
        nodeDegrees = adjMat.colwise().sum();
    }

    void SampleSomProblem::makeTedges(const SomUtils::MatD &T, SomUtils::MatD &Tedges, SomUtils::MatD &T1offset) const
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
        T1offset = T.col(src_);
    }

    void SampleSomProblem::makeTedges(const SomUtils::MatD &T, SomUtils::MatD &Tedges) const
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
        // ROFL_VAR2(T, Tedges)
    }

    void SampleSomProblem::POCRotateToMinimizeLastEntries(const SomUtils::MatD &x, SomUtils::MatD &Qtransp) const
    {
        Eigen::JacobiSVD<SomUtils::MatD> svd(x * x.transpose(), Eigen::ComputeFullV | Eigen::ComputeFullU); // TODO: ComputeFullV flag can probably be removed
        auto Q = svd.matrixU();
        Qtransp = Q.transpose();
    }

    void SampleSomProblem::makeTij1j2sEdges(int nodeId, const Eigen::ArrayXi &nodeDegrees, const SomUtils::MatD &Tedges,
                                            SomUtils::MatD &Tij1j2, SomUtils::MatD &Tij1j2tilde) const
    {
        ROFL_VAR1("makeTij1j2sEdges START")
        ROFL_VAR1(Tedges)

        // num_rows_T = size(T_edges, 1);
        int numRowsT = Tedges.rows();
        // % nrs = params.nrs;
        // d = params.d;
        auto nodeDeg = nodeDegrees(nodeId); // usually equal to low deg

        ROFL_ASSERT(Tij1j2.rows() == sz_.d_ && Tij1j2.cols() == nodeDeg)
        // SomUtils::MatD Tij1j2 = SomUtils::MatD::Zero(sz_.d_, nodeDeg);
        ROFL_ASSERT(Tij1j2tilde.rows() == numRowsT && Tij1j2tilde.cols() == nodeDeg)
        // SomUtils::MatD Tij1j2_tilde = SomUtils::MatD::Zero(num_rows_T, node_deg);

        Tij1j2.setZero();
        Tij1j2tilde.setZero();

        int found = 0;
        for (int e = 0; e < edges_.rows(); ++e)
        {
            auto eI = edges_(e, 0) - 1;
            // % e_j = edges(e, 2);
            if (eI == nodeId)
            {
                Tij1j2.col(found) = Tijs_.col(e);
                Tij1j2tilde.col(found) = -Tedges.col(e);
                found++;
            }
        }
        ROFL_ASSERT_VAR3(found == nodeDeg, nodeId, found, nodeDeg)

        ROFL_VAR1("makeTij1j2sEdges END")
        // ROFL_VAR2(Tij1j2, Tij1j2tilde)
    }

    void SampleSomProblem::orthComplement(const SomUtils::MatD &v, SomUtils::MatD &vOrth) const
    {
        // %Compute a basis for the orthogonal complement to the columns of v, i.e.,
        // %vOrth'*v=0

        // r=size(v,2);
        int r = v.cols();
        // [U,S,V]=svd(v);
        Eigen::JacobiSVD<SomUtils::MatD> svd(v, Eigen::ComputeFullV | Eigen::ComputeFullU); // TODO: ComputeFullV flag can probably be removed
        auto U = svd.matrixU();
        // vOrth=U(:,r+1:end);
        vOrth = U.block(0, r, U.rows(), U.cols() - r);
    }

    void SampleSomProblem::orthCompleteBasis(const SomUtils::MatD &Q, SomUtils::MatD &Qout) const
    {
        // %Complete the columns of Q (assumed orthonormal) to an orthonormal basis.
        Qout = Q;
        // [n,p]=size(Q);
        int n = Q.rows();
        int p = Q.cols();

        // [U,~,V]=svd([Q zeros(n,n-p)]);
        SomUtils::MatD Qpadded = SomUtils::MatD::Zero(n, n);
        Qpadded.block(0, 0, n, p) = Q;
        Eigen::JacobiSVD<SomUtils::MatD> svd(Qpadded, Eigen::ComputeFullV | Eigen::ComputeFullU); // TODO: ComputeFullV flag can probably be removed
        auto U = svd.matrixU();
        auto V = svd.matrixV();

        // Q=U*V';
        Qout = U * V.transpose();

        // %impose positive det(Q) (if Q is square) by flipping sign of last added
        // %column (if necessary)
        // if det(Q)<0 && p<n
        //     Q(:,end)=-Q(:,end);
        if (Qout.determinant() < 0 && p < n)
            Qout.col(Qout.cols() - 1) *= -1.0;
    }

    void SampleSomProblem::fliplr(const SomUtils::MatD &Ain, SomUtils::MatD &Aout) const
    {
        ROFL_ASSERT(Ain.rows() == Aout.rows() && Ain.cols() == Aout.cols())
        int ncols = Ain.cols();
        for (int i = 0; i < ncols; ++i)
        {
            Aout.col(i) = Ain.col(ncols - i - 1);
        }
    }

    void SampleSomProblem::flipud(const SomUtils::MatD &Ain, SomUtils::MatD &Aout) const
    {
        ROFL_ASSERT(Ain.rows() == Aout.rows() && Ain.cols() == Aout.cols())
        int nrows = Ain.rows();
        for (int i = 0; i < nrows; ++i)
        {
            Aout.row(i) = Ain.row(nrows - i - 1);
        }
    }

    void SampleSomProblem::align2d(const SomUtils::MatD &v, SomUtils::MatD &Qx) const
    {
        // Q=fliplr(orthComplement(v));
        SomUtils::MatD ocv = SomUtils::MatD::Zero(4, 2);
        orthComplement(v, ocv);
        SomUtils::MatD Q = SomUtils::MatD::Zero(ocv.rows(), ocv.cols());
        fliplr(ocv, Q);
        ROFL_VAR1("align2d first part")
        ROFL_VAR3(v, ocv, Q)

        // Qx=flipud(orthCompleteBasis(Q)');
        SomUtils::MatD ocb = SomUtils::MatD::Zero(4, 4);
        orthCompleteBasis(Q, ocb);
        Qx.resize(ocb.cols(), ocb.rows()); // NOT an error!!
        flipud(ocb.transpose(), Qx);
        ROFL_VAR1("align2d second part")
        ROFL_VAR2(ocb, Qx)
    }

    void SampleSomProblem::procrustesRb(const SomUtils::MatD &c, const SomUtils::MatD &q, SomUtils::MatD &RbEst) const
    {
        // [U,~,V]=svd(c*q');
        Eigen::JacobiSVD<SomUtils::MatD> svd(c * q.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV); // TODO: ComputeFullV flag can probably be removed
        auto U = svd.matrixU();
        auto V = svd.matrixV();

        // RbEst=U*diag([1 det(U*V')])*V';
        SomUtils::MatD tmp(SomUtils::MatD::Identity(2, 2));
        tmp(1, 1) = (U * V.transpose()).determinant();
        RbEst = U * tmp * V.transpose();

        ROFL_VAR1("procrustesRb output")
        ROFL_VAR3(c, q, RbEst);
    }

    void SampleSomProblem::recoverRiTilde(const SomUtils::MatD &RiTilde2, const SomUtils::MatD &Tijtilde,
                                          SomUtils::MatD &RiTildeEst1, SomUtils::MatD &RiTildeEst2) const
    {
        // PARAMS IN: Qx_edges*R_i_tilde2, Tij1j2_tilde
        ROFL_VAR1("Start of recoverRiTilde")
        ROFL_VAR2(RiTilde2, Tijtilde)

        // Qx = align2d(Tijtilde);
        SomUtils::MatD Qx(SomUtils::MatD::Zero(sz_.p_, sz_.p_));
        align2d(Tijtilde, Qx); // !! resizing of Qx happens inside align2d()

        // QxRitilde2Bot = Qx(3 : 4, :) * Ritilde2;
        ROFL_VAR2(Qx, RiTilde2)
        if (RiTilde2.rows() != Qx.cols() || Qx.rows() < 4)
        {
            ROFL_VAR1(".block() error avoided")
            RiTildeEst1 = RiTilde2;
            RiTildeEst2 = -RiTilde2;
            return;
        }
        auto QxRiTilde2bot = Qx.block(2, 0, 2, Qx.cols()) * RiTilde2;

        // [ U, ~, ~] = svd(QxRitilde2Bot, 'econ');
        Eigen::JacobiSVD<SomUtils::MatD> svd(QxRiTilde2bot, Eigen::ComputeThinU);
        auto U = svd.matrixU();

        // c = U( :, 2);
        auto c = U.col(1);

        // QLastRight = Qx(3 : 4, 4)';
        auto QLastRightT = Qx.block(2, 3, 2, 1);
        // auto QLastRight = QLastRightT.transpose();

        // RbEst = procrustesRb(c, QLastRight');
        SomUtils::MatD RbEst(SomUtils::MatD::Zero(2, 2));
        procrustesRb(c, QLastRightT, RbEst);

        // RitildeEst1 = Qx' * blkdiag(eye(2),-RbEst') * Qx * Ritilde2;
        SomUtils::MatD tmp1(SomUtils::MatD::Identity(4, 4));
        tmp1.block(2, 2, 2, 2) = -RbEst.transpose();
        RiTildeEst1 = Qx.transpose() * tmp1 * Qx * RiTilde2;

        // RitildeEst2=Qx'*blkdiag(eye(2),RbEst')*Qx*Ritilde2;
        SomUtils::MatD tmp2(SomUtils::MatD::Identity(4, 4));
        tmp2.block(2, 2, 2, 2) = RbEst.transpose();
        RiTildeEst2 = Qx.transpose() * tmp2 * Qx * RiTilde2;

        ROFL_VAR1("recoverRiTilde input")
        ROFL_VAR2(RiTilde2, Tijtilde)
        ROFL_VAR1("recoverRiTilde output")
        ROFL_VAR2(RiTildeEst1, RiTildeEst2)
    }

    void SampleSomProblem::dijkstraBT(int src, int n,
                                      const std::vector<int> &prev,
                                      std::vector<std::vector<int>> &list) const
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

    void SampleSomProblem::dijkstraBTedges(int src, int n,
                                           const std::vector<int> &prev, const Eigen::MatrixXi &edges,
                                           std::vector<std::vector<int>> &listEdges) const
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

    void SampleSomProblem::dijkstraSP(int n, int src,
                                      const Eigen::MatrixXi &adjmat,
                                      std::vector<double> &dist, std::vector<int> &prev) const
    {
        struct Vertex
        {
            int index;

            Vertex() : index(-1) {}

            Vertex(int idx) : index(idx) {}

            ~Vertex() {}
        };

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

    void SampleSomProblem::edgeDiffs2T(int src,
                                       const SomUtils::MatD &Tdiffs,
                                       int n,
                                       SomUtils::MatD &T) const
    {
        // nrs = size(T_diffs, 1);
        int nrs = Tdiffs.rows();
        ROFL_ASSERT(T.cols() == n)
        // booleans_T = boolean(0) * ones(N,1); % alg should stop when all these are 1
        Eigen::ArrayXi booleansT(Eigen::ArrayXi::Ones(n));
        // booleans_T(1) = boolean(1); % node 1 chosen as reference
        booleansT(0) = true;
        // T = zeros(nrs, N);
        T.setZero();
        // adjmat = edges2adjmatrix(edges);
        Eigen::MatrixXi adjMat(Eigen::MatrixXi::Zero(sz_.n_, sz_.n_));
        makeAdjMatFromEdges(adjMat);
        //
        int globalId = src;
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
                T.col(i) += Tdiffs.col(edgePath[j]);
            }
            T.col(i) *= -1;
            //     for ep = edge_path
            //         T(:,ii) = T(:,ii) + T_diffs(:,ep);
            //         disp('');
            //     end
            //     T(:,ii) = -T(:,ii);
        }
    }

    bool SampleSomProblem::recoverySEdN(int staircaseStepIdx,
                                        const SomUtils::VecMatD &RmanoptOut, const SomUtils::MatD &TmanoptOut,
                                        SomUtils::VecMatD &Rrecovered, SomUtils::MatD &Trecovered)
    {
        rsRecoverySuccess_ = true;
        ROFL_ASSERT(Rrecovered.size() == sz_.n_ && Trecovered.rows() == sz_.d_ && Trecovered.cols() == sz_.n_)

        // !! In the way recovery on SE(d)^N is formulated now, we have to perform it even if nrs = d

        int nrs = staircaseStepIdx - 1;
        nrs = TmanoptOut.rows();
        int lowDeg = 2; // TODO: not necessarily 2 in more complex graph cases (?)

        Eigen::ArrayXi nodeDegrees(Eigen::ArrayXi::Zero(sz_.n_));
        computeNodeDegrees(nodeDegrees);

        // auto nodesHighDeg = nodeDegrees > lowDeg;
        Eigen::ArrayXi nodesHighDeg(Eigen::ArrayXi::Zero(sz_.n_));
        Eigen::ArrayXi nodesLowDeg(Eigen::ArrayXi::Zero(sz_.n_));

        for (int i = 0; i < sz_.n_; ++i)
        {
            if (nodeDegrees(i, 0) > lowDeg)
                nodesHighDeg(i, 0) = 1;
            else
                nodesLowDeg(i, 0) = 1;
        }
        ROFL_VAR1(nodeDegrees.transpose());
        ROFL_VAR1(nodesHighDeg.transpose());
        ROFL_VAR1(nodesLowDeg.transpose())

        // [T_edges, T1_offset] = make_T_edges(T_manopt_out, edges);
        SomUtils::MatD Tedges(SomUtils::MatD::Zero(nrs, numEdges_));
        makeTedges(TmanoptOut, Tedges);

        // RT_stacked_high_deg = [ matStackH(R_manopt_out( :, :, nodes_high_deg)), T_edges ];
        auto numNodesHighDeg = nodesHighDeg.sum();
        ROFL_VAR1(numNodesHighDeg);
        SomUtils::MatD RTstackedHighDeg(SomUtils::MatD::Zero(nrs, sz_.d_ * numNodesHighDeg + Tedges.cols()));
        SomUtils::MatD RstackedHighDeg(SomUtils::MatD::Zero(nrs, sz_.d_ * numNodesHighDeg));
        SomUtils::VecMatD RmanoptOutHighDeg;
        for (int i = 0; i < sz_.n_; ++i)
        {
            if (nodesHighDeg(i, 0) != 0)
                RmanoptOutHighDeg.push_back(RmanoptOut[i]);
        }
        ROFL_VAR1("hstack call from here");
        SomUtils::hstack(RmanoptOutHighDeg, RstackedHighDeg);
        RTstackedHighDeg.block(0, 0, nrs, sz_.d_ * numNodesHighDeg) = RstackedHighDeg;
        RTstackedHighDeg.block(0, sz_.d_ * numNodesHighDeg, nrs, numEdges_) = Tedges;

        // ROFL_VAR1(RTstackedHighDeg)

        SomUtils::MatD QxEdges(SomUtils::MatD::Zero(nrs, nrs));
        POCRotateToMinimizeLastEntries(RTstackedHighDeg, QxEdges);

        // ROFL_VAR1(QxEdges)

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
        // for (int i = 0; i < numNodesHighDeg; ++i)
        //     ROFL_VAR1(Rtilde2edges[i])

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

        // for (int i = 0; i < sz_.n_; ++i)
        //     ROFL_VAR1(Rrecovered[i])

        if (!nodesLowDeg.any())
        {
            ROFL_VAR1("No nodes low deg!");
            SomUtils::MatD TdiffsShifted = QxEdges * Tedges; // this has last row to 0
            edgeDiffs2T(src_, TdiffsShifted.block(0, 0, sz_.d_, TdiffsShifted.cols()), sz_.n_, Trecovered);
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
                    ROFL_VAR1("hstack call from here");
                    SomUtils::hstack(Rgt_, RgtSt);
                    ROFL_VAR1(RgtSt);

                    Xgt.block(0, 0, sz_.d_, RgtSt.cols()) = RgtSt;
                    Xgt.block(0, RgtSt.cols(), sz_.d_, Tgt_.cols()) = Tgt_;

                    // disp("cost_gt")
                    // disp(cost_gt)
                    double costGt = costEigen(Rgt_, Tgt_);
                    ROFL_VAR1(costGt)

                    SomUtils::MatD XmanoptOut(SomUtils::MatD::Zero(nrs, sz_.d_ * sz_.n_ + sz_.d_ * sz_.n_));
                    SomUtils::MatD RmanoptOutSt(SomUtils::MatD::Zero(nrs, sz_.d_ * sz_.n_));
                    ROFL_VAR1("hstack call from here");
                    SomUtils::hstack(RmanoptOut, RmanoptOutSt);
                    ROFL_VAR1(RmanoptOutSt);

                    XmanoptOut.block(0, 0, nrs, RmanoptOutSt.cols()) = RmanoptOutSt;
                    XmanoptOut.block(0, RmanoptOutSt.cols(), nrs, TmanoptOut.cols()) = TmanoptOut;

                    double costManoptOutput = costEigen(RmanoptOut, TmanoptOut);

                    // disp("cost_manopt_output")
                    // disp(cost_manopt_output)
                    ROFL_VAR1(costManoptOutput);
                    // T_diffs_shifted = Qx_edges * T_edges;
                    auto TdiffsShifted = QxEdges * Tedges; // this has last row to 0
                    ROFL_VAR3(TdiffsShifted, QxEdges, Tedges)

                    // [~, Tij1j2_tilde] = make_Tij1j2s_edges(node_id, T_diffs_shifted, Tijs, edges, params);
                    SomUtils::MatD Tij1j2(SomUtils::MatD::Zero(sz_.d_, nodeDeg));
                    SomUtils::MatD Tij1j2tilde(SomUtils::MatD::Zero(nrs, nodeDeg));
                    makeTij1j2sEdges(nodeId, nodeDegrees, TdiffsShifted, Tij1j2, Tij1j2tilde);
                    // ROFL_VAR2(Tij1j2, Tij1j2tilde)

                    // [ RitildeEst1, RitildeEst2, ~, ~] = recoverRitilde(Qx_edges * R_i_tilde2, Tij1j2_tilde);
                    SomUtils::MatD RiTildeEst1(SomUtils::MatD::Zero(nrs, sz_.d_));
                    SomUtils::MatD RiTildeEst2(SomUtils::MatD::Zero(nrs, sz_.d_));
                    // ROFL_VAR2(QxEdges, RiTilde2)
                    recoverRiTilde(QxEdges * RiTilde2, Tij1j2tilde, RiTildeEst1, RiTildeEst2); // TODO: add possibility of returning "local" Qx
                    // ROFL_VAR2(Tij1j2, Tij1j2tilde)
                    // ROFL_VAR2(RiTildeEst1, RiTildeEst2)

                    // disp('')
                    std::cout << std::endl; // TODO : how to decide between RitildeEst1, RitildeEst2 ? ? det_RitildeEst1 = det(RitildeEst1(1 : d, :));
                    // det_RitildeEst2 = det(RitildeEst2(1 : d, :));
                    auto detRiTildeEst1 = RiTildeEst1.block(0, 0, sz_.d_, sz_.d_).determinant();
                    auto detRiTildeEst2 = RiTildeEst2.block(0, 0, sz_.d_, sz_.d_).determinant();
                    ROFL_VAR2(detRiTildeEst1, detRiTildeEst2)

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
                        rsRecoverySuccess_ = false; // maybe add possibility to return this also
                        // save('data/zerodet_ws.mat')
                        // ROFL_ASSERT(0) // TODO: add this line for non-noisy cases after recoverRiTilde is implemented
                        // }
                    }
                    // T_recovered = edge_diffs_2_T(T_diffs_shifted(1 : d, :), edges, N);
                    edgeDiffs2T(src_, TdiffsShifted.block(0, 0, sz_.d_, TdiffsShifted.cols()), sz_.n_, Trecovered);
                    std::cout << std::endl;
                }
            }
        }

        // checking that cost has not changed during "recovery" X_recovered.T = T_recovered;
        // X_recovered.R = R_recovered;
        // cost_out = rsom_cost_base(X_recovered, problem_struct_next);
        // disp("cost_out")
        // disp(cost_out)
        SomUtils::MatD Xrecovered(SomUtils::MatD::Zero(sz_.d_, sz_.d_ * sz_.n_ + sz_.d_ * sz_.n_));
        SomUtils::MatD RrecoveredSt(SomUtils::MatD::Zero(sz_.d_, sz_.d_ * sz_.n_));
        ROFL_VAR1("hstack call from here");
        SomUtils::hstack(Rrecovered, RrecoveredSt);
        ROFL_VAR1(RrecoveredSt);
        Xrecovered.block(0, 0, sz_.d_, RrecoveredSt.cols()) = RrecoveredSt;
        Xrecovered.block(0, RrecoveredSt.cols(), sz_.d_, Trecovered.cols()) = Trecovered;

        ROFL_VAR1("Before end of recoverySEdN()")
        //  disp("[matStackH(X_gt.R); matStackH(R_recovered)]");
        //  disp([matStackH(X_gt.R); matStackH(R_recovered)]);
        for (int i = 0; i < sz_.n_; ++i)
        {
            ROFL_VAR3(i, Rgt_[i], Rrecovered[i])
            ROFL_VAR2(i, Tgt_.col(i).transpose())
            ROFL_VAR2(TmanoptOut.col(i).transpose(), Trecovered.col(i).transpose())
        }
        return rsRecoverySuccess_;
    }

    bool SampleSomProblem::globalize(int src, const SomUtils::VecMatD &Rsedn, const SomUtils::MatD &Tsedn,
                                     SomUtils::VecMatD &Rout, SomUtils::MatD &Tout)
    {
        // R_recovered -> Rsedn
        // T_recovered -> Tsedn

        // GLOBALIZATION! -> probably put it in another function?
        // R_global = R_recovered(:,:,1) * X_gt.R(:,:,1)'; %!!
        auto Rglobal = Rsedn[src] * Rgt_[src].transpose();
        ROFL_VAR1(Rglobal)

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
            ROFL_VAR4(i, Rgt_[i], RrecoveredGlobal[i], SomUtils::isEqualFloats(Rgt_[i], RrecoveredGlobal[i]))
        }

        // T_global = R_global * T_recovered(:,1) - X_gt.T(:,1); %!!
        auto Tglobal = Rglobal * Tsedn.col(src) - Tgt_.col(src);
        ROFL_VAR4(Rglobal, (Rglobal * Tsedn.col(src)).transpose(), Tsedn.col(src).transpose(), Tgt_.col(src).transpose());
        ROFL_VAR2(Tglobal.transpose(), Tsedn.col(src).transpose());
        // code for making all translation global at once
        // disp("[X_gt.T; T_recovered]");

        // T_recovered_global = R_global' * T_recovered - T_global;
        // disp([X_gt.T; T_recovered_global]);
        ROFL_VAR3(Rglobal.transpose(), Tsedn, Tglobal)
        SomUtils::MatD TglobalRepmat(SomUtils::MatD::Zero(sz_.d_, sz_.n_));
        for (int i = 0; i < sz_.n_; ++i)
        {
            TglobalRepmat.col(i) = Tglobal; // TODO: maybe use some other adv init
        }
        auto TrecoveredGlobal = Rglobal.transpose() * Tsedn - TglobalRepmat;
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
            if (!SomUtils::isEqualFloats(RgtI, RrecovIglobal))
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
            ROFL_VAR4(i, TgtI.transpose(), Tsedn.col(i).transpose(), TrecovIglobal.transpose());

            //     disp("is_equal_floats(T_gt_i, T_recov_i_global)")
            //     disp(is_equal_floats(T_gt_i, T_recov_i_global))
            ROFL_VAR1(SomUtils::isEqualFloats(TgtI, TrecovIglobal))
            //     if (~is_equal_floats(T_gt_i, T_recov_i_global))
            // %         error("transl found NOT equal")
            //         fprintf("ERROR in recovery: T_GLOBAL\n");
            if (!SomUtils::isEqualFloats(TgtI, TrecovIglobal, 1e-3))
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
        ROFL_VAR1("hstack call from here");
        SomUtils::hstack(RrecoveredGlobal, RoutSt);
        ROFL_VAR1(RoutSt);
        Xout.block(0, 0, sz_.d_, RoutSt.cols()) = RoutSt;
        Xout.block(0, RoutSt.cols(), sz_.d_, Tout.cols()) = TrecoveredGlobal;

        Rout = RrecoveredGlobal;
        Tout = TrecoveredGlobal;
        ROFL_VAR1(costEigen(Rout, Tout))

        // transf_out = RT2G(R_recovered_global, T_recovered_global); %rsom_genproc() function output

        // DETERMINANTS CHECK
        std::vector<double> multidetRrecovered, multidetRrecoveredGlobal;
        // disp('multidet(R_recovered)')
        // disp(multidet(R_recovered))
        SomUtils::multidet(Rsedn, multidetRrecovered);
        for (int i = 0; i < sz_.n_; ++i)
            ROFL_VAR2(i, multidetRrecovered[i]);

        // disp('multidet(R_recovered_global)')
        // disp(multidet(R_recovered_global))
        SomUtils::multidet(Rout, multidetRrecoveredGlobal);
        for (int i = 0; i < sz_.n_; ++i)
            ROFL_VAR2(i, multidetRrecoveredGlobal[i]);

        return rsRecoverySuccess_;
    }

} // end of namespace ROPTLIB