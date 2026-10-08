#include "problems_Lsom.h"

namespace ROPTLIB
{
    LsomProblem::LsomProblem() {}

    LsomProblem::LsomProblem(SomUtils::SomSize somSz, SomUtils::MatD &tijs, Eigen::MatrixXi &edges)
    {
        sz_ = somSz;
        tijs_ = tijs;
        edges_ = edges;
        numEdges_ = tijs_.cols();
        fullSz_ = sz_.d_ * sz_.p_ * sz_.n_ + sz_.p_ * sz_.n_ + numEdges_;

        Rgt_.resize(sz_.n_, SomUtils::MatD::Zero(sz_.d_, sz_.d_));
        Tgt_ = SomUtils::MatD::Zero(sz_.d_, sz_.n_);
        LambdasGt_ = SomUtils::MatD::Zero(numEdges_, 1);

        rho_ = 0.0; // TODO: add rho_ as input parameter in another constructor

        src_ = 0; // TODO: src_ VS src (for sure in globalize, maybe also in other places)

        costCurr_ = 1e+10;

        usePIM_ = true; // default to PIM; can be changed with setter if needed

        pimMaxIterations_ = 5000; // default max iterations for PIM; can be changed with setter if needed

        reluScaleCompensation_ = false; // default to false; can be changed with setter if needed

        a_ = 2.0; // default value for log-based compensation; can be changed with setter if needed

        b_ = -a_ / ((a_ - 1) * (a_ - 1)); // dependent on a_

        enableRs_ = true; // default to false; can be changed with setter if needed

        maxIterAdmm_ = 10000;  // default value, can be changed by setter
        tolAdmmPrimal_ = 1e-8; // default value, can be changed by setter
        tolAdmmDual_ = 1e-8;   // default value, can be changed by setter

        // zAdmm_ = SomUtils::MatD::Ones(numEdges_, 1);
        if (reluScaleCompensation_)
            zAdmm_ = 5.0 * SomUtils::MatD::Ones(numEdges_, 1);
        else
            zAdmm_ = 10.0 * SomUtils::MatD::Ones(numEdges_, 1);

        yAdmm_ = SomUtils::MatD::Zero(numEdges_, 1);
        muAdmm_ = 0.5; // default value, can be changed by setter

        ssomInitguess_ = true;

        firstZadmmLambdas_ = false;

        performGlobalization_ = true; // default to false; can be changed with setter if needed
    }

    LsomProblem::LsomProblem(const SomUtils::SomSize somSz, const SomUtils::MatD &tijs, const Eigen::MatrixXi &edges)
    {
        sz_ = somSz;
        tijs_ = tijs;
        edges_ = edges;
        numEdges_ = tijs_.cols();
        fullSz_ = sz_.d_ * sz_.p_ * sz_.n_ + sz_.p_ * sz_.n_ + numEdges_;

        Rgt_.resize(sz_.n_, SomUtils::MatD::Zero(sz_.d_, sz_.d_));
        Tgt_ = SomUtils::MatD::Zero(sz_.d_, sz_.n_);
        LambdasGt_ = SomUtils::MatD::Zero(numEdges_, 1);

        rho_ = 0.0; // TODO: add rho_ as input parameter in another constructor

        src_ = 0;

        costCurr_ = 1e+10;

        usePIM_ = true; // default to PIM; can be changed with a setter if needed

        pimMaxIterations_ = 5000; // default max iterations for PIM; can be changed with a setter if needed

        reluScaleCompensation_ = false; // default to false; can be changed with setter if needed

        a_ = 2.0; // default value for log-based compensation; can be changed with setter if needed

        b_ = -a_ / ((a_ - 1) * (a_ - 1)); // dependent on a_

        enableRs_ = true; // default to false; can be changed with setter if needed

        maxIterAdmm_ = 10000;  // default value, can be changed by setter
        tolAdmmPrimal_ = 1e-8; // default value, can be changed by setter
        tolAdmmDual_ = 1e-8;   // default value, can be changed by setter

        // zAdmm_ = SomUtils::MatD::Ones(numEdges_, 1);
        if (reluScaleCompensation_)
            zAdmm_ = 5.0 * SomUtils::MatD::Ones(numEdges_, 1);
        else
            zAdmm_ = 10.0 * SomUtils::MatD::Ones(numEdges_, 1);

        yAdmm_ = SomUtils::MatD::Zero(numEdges_, 1);
        muAdmm_ = 0.5; // default value, can be changed by setter

        ssomInitguess_ = true;

        firstZadmmLambdas_ = false;

        performGlobalization_ = true; // default to false; can be changed with setter if needed
    }

    LsomProblem::~LsomProblem() {};

    realdp LsomProblem::f(const Variable &x) const
    {
        SomUtils::MatD xEigen(fullSz_, 1);
        RoptToEig(x, xEigen);
        // ROFL_VAR1(x);

        realdp cost = 1e+10;
        if (reluScaleCompensation_)
            cost = costEigenVecRelu(xEigen);
        else
            cost = costEigenVec(xEigen);

        ROFL_VAR2("cost in LsomProblem::f()", cost);
        // ROFL_ASSERT(!std::isnan(corr));
        if (std::isnan(cost))
        {
            ROFL_ERR("Cost is NaN")
            ROFL_VAR1(reluScaleCompensation_)

            double costDebug = 0.0f;

            for (int e = 0; e < numEdges_; ++e)
            {
                SomUtils::MatD Ri(SomUtils::MatD::Zero(sz_.p_, sz_.d_));
                SomUtils::MatD Ti(SomUtils::MatD::Zero(sz_.p_, 1));
                SomUtils::MatD Tj(SomUtils::MatD::Zero(sz_.p_, 1));
                double lambdaE = 0.0;

                SomUtils::VecD tij(SomUtils::VecD::Zero(sz_.d_));
                tij = tijs_.col(e);

                int i = edges_(e, 0) - 1; // !! -1
                int j = edges_(e, 1) - 1; // !! -1
                getRi(xEigen, Ri, i);
                getTi(xEigen, Ti, i);
                getTi(xEigen, Tj, j);
                getLambdaI(xEigen, lambdaE, e);

                // ROFL_VAR3(i, j, e);
                // ROFL_VAR4(Ri, tij.transpose(), Ti.transpose(), Tj.transpose());

                auto a = Ti - Tj;
                auto b = Ri * tij;
                auto costLambdaEe = (a.transpose() * a + 2 * lambdaE * (a.transpose() * b) + lambdaE * lambdaE * (b.transpose() * b)).trace();

                // l = lambda_e;
                // if l<1
                //     scale_compensation_ee=-1/a_log*log(a_log*l-1)...
                //         +1/(a_log-1)*(l-1)...
                //         +b_log/2*(l-1)^2;
                // else
                //     scale_compensation_ee=0;

                double scaleCompensation = 0.0;
                if (lambdaE <= 1 / a_)
                {
                    if (fabs(rho_) < 1e-6)
                    {
                        scaleCompensation = 0;
                    }
                    else
                    {
                        scaleCompensation = std::nan("nan");
                    }
                }
                else if (lambdaE < 1.0)
                {
                    scaleCompensation = (-1.0 / a_) * log(a_ * lambdaE - 1.0) + (1.0 / (a_ - 1.0)) * (lambdaE - 1.0) + (b_ / 2.0) * (lambdaE - 1.0) * (lambdaE - 1.0);
                }

                ROFL_VAR5(e, scaleCompensation, log(a_ * lambdaE - 1.0), a_, lambdaE)

                if (std::isnan(lambdaE))
                {
                    ROFL_ERR("lambdaE is NaN")
                    ROFL_VAR5(e, a_, b_, lambdaE, costLambdaEe)
                    ROFL_VAR3(a.transpose(), b.transpose(), xEigen.transpose())
                    ROFL_ASSERT(0)
                }

                // ROFL_VAR6(e, a_, b_, lambdaE, costLambdaEe, scaleCompensation);

                costDebug += costLambdaEe + rho_ * scaleCompensation;
            }
            ROFL_VAR1(costDebug)

            ROFL_ASSERT(0)
        }

        ROFL_ASSERT(!std::isnan(cost));

        // Vector *resultEgrad;
        // *resultEgrad = Domain->RandominManifold();
        // EucGrad(x, resultEgrad);
        // x.AddToFields("EGrad", *resultEgrad); //x should have been const?? maybe only its reference

        return cost; // checked -> the - here should be OK
    };

    double LsomProblem::costEigenVecReluSEdN(const SomUtils::MatD &xEigen) const
    {
        double cost = 0.0f;
        for (int e = 0; e < numEdges_; ++e)
        {
            SomUtils::MatD Ri(SomUtils::MatD::Zero(sz_.d_, sz_.d_));
            SomUtils::MatD Ti(SomUtils::MatD::Zero(sz_.d_, 1));
            SomUtils::MatD Tj(SomUtils::MatD::Zero(sz_.d_, 1));

            SomUtils::VecD tij(SomUtils::VecD::Zero(sz_.d_));
            tij = tijs_.col(e);

            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1
            getRiSEdN(xEigen, Ri, i);
            getTiSEdN(xEigen, Ti, i);
            getTiSEdN(xEigen, Tj, j);
            double lambdaE = 0.0;
            getLambdaISEdN(xEigen, lambdaE, e);

            // ROFL_VAR3(i, j, e);
            // ROFL_VAR4(Ri, tij.transpose(), Ti.transpose(), Tj.transpose());

            auto a = Ti - Tj;
            auto b = Ri * tij;
            auto costLambdaEe = (a.transpose() * a + 2 * lambdaE * (a.transpose() * b) + lambdaE * lambdaE * (b.transpose() * b)).trace();
            auto costReluEe = SomUtils::ReLU(lsomReLUargument(lambdaE));
            cost += costLambdaEe + rho_ * costReluEe * costReluEe;
        }

        // cost_out = cost_out + y'*(vec(z)-vec(lambdas))+0.5 * mu * norm(vec(z)-vec(lambdas))^2;
        SomUtils::MatD LambdasEigen(numEdges_, 1);
        getScales(xEigen, LambdasEigen);
        ROFL_VAR2(zAdmm_.transpose(), LambdasEigen.transpose());
        auto admmCost1 = yAdmm_.transpose() * (zAdmm_ - LambdasEigen);
        ROFL_ASSERT_VAR2(admmCost1.rows() == 1 && admmCost1.cols() == 1, admmCost1.rows(), admmCost1.cols());
        double admmCost = admmCost1(0, 0) + 0.5 * muAdmm_ * ((zAdmm_ - LambdasEigen).squaredNorm());
        return cost + admmCost;
    }

    double LsomProblem::costEigenVecSEdN(const SomUtils::MatD &xEigen) const
    {
        double cost = 0.0f;
        for (int e = 0; e < numEdges_; ++e)
        {
            SomUtils::MatD Ri(SomUtils::MatD::Zero(sz_.d_, sz_.d_));
            SomUtils::MatD Ti(SomUtils::MatD::Zero(sz_.d_, 1));
            SomUtils::MatD Tj(SomUtils::MatD::Zero(sz_.d_, 1));

            SomUtils::VecD tij(SomUtils::VecD::Zero(sz_.d_));
            tij = tijs_.col(e);

            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1
            getRiSEdN(xEigen, Ri, i);
            getTiSEdN(xEigen, Ti, i);
            getTiSEdN(xEigen, Tj, j);
            double lambdaE = 0.0;
            getLambdaISEdN(xEigen, lambdaE, e);

            // ROFL_VAR3(i, j, e);
            // ROFL_VAR4(Ri, tij.transpose(), Ti.transpose(), Tj.transpose());

            auto a = Ti - Tj;
            auto b = Ri * tij;
            auto costLambdaEe = (a.transpose() * a + 2 * lambdaE * (a.transpose() * b) + lambdaE * lambdaE * (b.transpose() * b)).trace();

            // l = lambda_e;
            // if l<1
            //     scale_compensation_ee=-1/a_log*log(a_log*l-1)...
            //         +1/(a_log-1)*(l-1)...
            //         +b_log/2*(l-1)^2;
            // else
            //     scale_compensation_ee=0;

            double scaleCompensation = 0.0;

            if (lambdaE <= 1 / a_)
            {
                // scaleCompensation = std::nan("nan");
                scaleCompensation = 1e+10;
            }
            else if (lambdaE < 1.0)
            {
                scaleCompensation = (-1.0 / a_) * log(a_ * lambdaE - 1.0) + (1.0 / (a_ - 1.0)) * (lambdaE - 1.0) + (b_ / 2.0) * (lambdaE - 1.0) * (lambdaE - 1.0);
            }

            cost += costLambdaEe + rho_ * scaleCompensation;
        }
        SomUtils::MatD LambdasEigen(numEdges_, 1);
        getScales(xEigen, LambdasEigen);
        auto admmCost1 = yAdmm_.transpose() * (zAdmm_ - LambdasEigen);
        ROFL_ASSERT_VAR2(admmCost1.rows() == 1 && admmCost1.cols() == 1, admmCost1.rows(), admmCost1.cols());
        double admmCost = admmCost1(0, 0) + 0.5 * muAdmm_ * ((zAdmm_ - LambdasEigen).squaredNorm());
        return cost + admmCost;
    }

    double LsomProblem::costEigenRelu(const SomUtils::VecMatD &Reigen, const SomUtils::MatD &Teigen, const SomUtils::MatD &LambdasEigen) const
    {
        double cost = 0.0f;
        for (int e = 0; e < numEdges_; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1

            auto Ri = Reigen[i];
            auto Ti = Teigen.col(i);
            auto Tj = Teigen.col(j);
            auto lambdaE = LambdasEigen(e, 0);

            SomUtils::VecD tij(SomUtils::VecD::Zero(sz_.d_));
            tij = tijs_.col(e);

            // ROFL_VAR3(i, j, e);
            // ROFL_VAR5(Ri, tij.transpose(), Ti.transpose(), Tj.transpose(), lambdaE);

            auto a = Ti - Tj;
            auto b = Ri * tij;
            auto costLambdaEe = (a.transpose() * a + 2 * lambdaE * (a.transpose() * b) + lambdaE * lambdaE * (b.transpose() * b)).trace();
            auto costReluEe = SomUtils::ReLU(lsomReLUargument(lambdaE));
            cost += costLambdaEe + rho_ * costReluEe * costReluEe;

            // ROFL_VAR1(lambdaE);
            // ROFL_VAR5(e, a.transpose(), b.transpose(), costLambdaEe, costReluEe);
        }

        auto admmCost1 = yAdmm_.transpose() * (zAdmm_ - LambdasEigen);
        ROFL_ASSERT_VAR2(admmCost1.rows() == 1 && admmCost1.cols() == 1, admmCost1.rows(), admmCost1.cols());
        double admmCost = admmCost1(0, 0) + 0.5 * muAdmm_ * ((zAdmm_ - LambdasEigen).squaredNorm());
        return cost + admmCost;
    }

    double LsomProblem::costEigen(const SomUtils::VecMatD &Reigen, const SomUtils::MatD &Teigen, const SomUtils::MatD &LambdasEigen) const
    {
        double cost = 0.0f;
        for (int e = 0; e < numEdges_; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1

            auto Ri = Reigen[i];
            auto Ti = Teigen.col(i);
            auto Tj = Teigen.col(j);
            auto lambdaE = LambdasEigen(e, 0);

            SomUtils::VecD tij(SomUtils::VecD::Zero(sz_.d_));
            tij = tijs_.col(e);

            // ROFL_VAR3(i, j, e);
            // ROFL_VAR5(Ri, tij.transpose(), Ti.transpose(), Tj.transpose(), lambdaE);

            auto a = Ti - Tj;
            auto b = Ri * tij;
            auto costLambdaEe = (a.transpose() * a + 2 * lambdaE * (a.transpose() * b) + lambdaE * lambdaE * (b.transpose() * b)).trace();

            // l = lambda_e;
            // if l<1
            //     scale_compensation_ee=-1/a_log*log(a_log*l-1)...
            //         +1/(a_log-1)*(l-1)...
            //         +b_log/2*(l-1)^2;
            // else
            //     scale_compensation_ee=0;

            double scaleCompensation = 0.0;
            if (lambdaE <= 1 / a_)
            {
                scaleCompensation = 1e+10;
            }
            else if (lambdaE < 1.0)
            {
                scaleCompensation = (-1.0 / a_) * log(a_ * lambdaE - 1.0) + (1.0 / (a_ - 1.0)) * (lambdaE - 1.0) + (b_ / 2.0) * (lambdaE - 1.0) * (lambdaE - 1.0);
            }

            // ROFL_VAR1(lambdaE);
            // ROFL_VAR5(e, a.transpose(), b.transpose(), costLambdaEe, scaleCompensation);

            cost += costLambdaEe + rho_ * scaleCompensation;
        }
        // cost_out = cost_out + y'*(vec(z)-vec(lambdas))+0.5 * mu * norm(vec(z)-vec(lambdas))^2
        auto admmCost1 = yAdmm_.transpose() * (zAdmm_ - LambdasEigen);
        ROFL_ASSERT_VAR2(admmCost1.rows() == 1 && admmCost1.cols() == 1, admmCost1.rows(), admmCost1.cols());
        double admmCost = admmCost1(0, 0) + 0.5 * muAdmm_ * ((zAdmm_ - LambdasEigen).squaredNorm());
        return cost + admmCost;
    }

    double LsomProblem::costEigenVecRelu(const SomUtils::MatD &xEigen) const
    {
        double cost = 0.0f;
        for (int e = 0; e < numEdges_; ++e)
        {
            SomUtils::MatD Ri(SomUtils::MatD::Zero(sz_.p_, sz_.d_));
            SomUtils::MatD Ti(SomUtils::MatD::Zero(sz_.p_, 1));
            SomUtils::MatD Tj(SomUtils::MatD::Zero(sz_.p_, 1));
            double lambdaE = 0.0;

            SomUtils::VecD tij(SomUtils::VecD::Zero(sz_.d_));
            tij = tijs_.col(e);

            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1
            getRi(xEigen, Ri, i);
            getTi(xEigen, Ti, i);
            getTi(xEigen, Tj, j);
            getLambdaI(xEigen, lambdaE, e);

            // ROFL_VAR3(i, j, e);
            // ROFL_VAR4(Ri, tij.transpose(), Ti.transpose(), Tj.transpose());

            auto a = Ti - Tj;
            auto b = Ri * tij;
            auto costLambdaEe = (a.transpose() * a + 2 * lambdaE * (a.transpose() * b) + lambdaE * lambdaE * (b.transpose() * b)).trace();
            auto costReluEe = SomUtils::ReLU(lsomReLUargument(lambdaE));
            cost += costLambdaEe + rho_ * costReluEe * costReluEe;
        }
        SomUtils::MatD LambdasEigen(numEdges_, 1);
        getScales(xEigen, LambdasEigen); // TODO: can be called before for loop to save some time, but should be fine for now
        // cost_out = cost_out + y'*(vec(z)-vec(lambdas))+0.5 * mu * norm(vec(z)-vec(lambdas))^2
        ROFL_VAR1(yAdmm_.transpose());
        ROFL_VAR1(zAdmm_.transpose())
        ROFL_VAR1(LambdasEigen.transpose());

        auto admmCost1 = yAdmm_.transpose() * (zAdmm_ - LambdasEigen);
        ROFL_ASSERT_VAR2(admmCost1.rows() == 1 && admmCost1.cols() == 1, admmCost1.rows(), admmCost1.cols());
        ROFL_VAR1(admmCost1);
        ROFL_ASSERT_VAR1(admmCost1(0, 0) >= 0.0, admmCost1);
        double admmCost = admmCost1(0, 0) + 0.5 * muAdmm_ * ((zAdmm_ - LambdasEigen).squaredNorm());
        ROFL_VAR2(cost, admmCost);
        ROFL_ASSERT_VAR1(admmCost >= 0.0, admmCost);

        return cost + admmCost;
    }

    double LsomProblem::costEigenVec(const SomUtils::MatD &xEigen) const
    {
        double cost = 0.0f;

        for (int e = 0; e < numEdges_; ++e)
        {
            SomUtils::MatD Ri(SomUtils::MatD::Zero(sz_.p_, sz_.d_));
            SomUtils::MatD Ti(SomUtils::MatD::Zero(sz_.p_, 1));
            SomUtils::MatD Tj(SomUtils::MatD::Zero(sz_.p_, 1));
            double lambdaE = 0.0;

            SomUtils::VecD tij(SomUtils::VecD::Zero(sz_.d_));
            tij = tijs_.col(e);

            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1
            getRi(xEigen, Ri, i);
            getTi(xEigen, Ti, i);
            getTi(xEigen, Tj, j);
            getLambdaI(xEigen, lambdaE, e);

            // ROFL_VAR3(i, j, e);
            // ROFL_VAR4(Ri, tij.transpose(), Ti.transpose(), Tj.transpose());

            auto a = Ti - Tj;
            auto b = Ri * tij;
            auto costLambdaEe = (a.transpose() * a + 2 * lambdaE * (a.transpose() * b) + lambdaE * lambdaE * (b.transpose() * b)).trace();

            // l = lambda_e;
            // if l<1
            //     scale_compensation_ee=-1/a_log*log(a_log*l-1)...
            //         +1/(a_log-1)*(l-1)...
            //         +b_log/2*(l-1)^2;
            // else
            //     scale_compensation_ee=0;

            double scaleCompensation = 0.0;
            if (lambdaE <= 1 / a_)
            {
                scaleCompensation = 1e+10;
            }
            else if (lambdaE < 1.0)
            {
                scaleCompensation = (-1.0 / a_) * log(a_ * lambdaE - 1.0) + (1.0 / (a_ - 1.0)) * (lambdaE - 1.0) + (b_ / 2.0) * (lambdaE - 1.0) * (lambdaE - 1.0);
            }

            // ROFL_VAR6(e, a_, b_, lambdaE, costLambdaEe, scaleCompensation);

            cost += costLambdaEe + rho_ * scaleCompensation;
        }
        SomUtils::MatD LambdasEigen(numEdges_, 1);
        getScales(xEigen, LambdasEigen); // TODO: can be called before for loop to save some time, but should be fine for now
        // cost_out = cost_out + y'*(vec(z)-vec(lambdas))+0.5 * mu * norm(vec(z)-vec(lambdas))^2
        auto admmCost1 = yAdmm_.transpose() * (zAdmm_ - LambdasEigen);
        ROFL_ASSERT_VAR2(admmCost1.rows() == 1 && admmCost1.cols() == 1, admmCost1.rows(), admmCost1.cols());
        double admmCost = admmCost1(0, 0) + 0.5 * muAdmm_ * ((zAdmm_ - LambdasEigen).squaredNorm());
        ROFL_VAR2(cost, admmCost);
        return cost + admmCost;
    }

    // Vector &LsomProblem::EucGrad(const Variable &x, Vector *result) const
    // {
    // };

    Vector &LsomProblem::RieGrad(const Variable &x, Vector *result) const
    {
        // result->NewMemoryOnWrite();
        // result->GetElement(0) = x.Field("B1x1D1");
        // result->GetElement(1) = x.Field("B2x2D2");
        // result->GetElement(2) = x.Field("B3x3D3");
        // Domain->ScalarTimesVector(x, 2, *result, result);

        // result->Print("RieGrad: printing result at start of function (should be empty)");

        // ROFL_VAR1(sz_.p_)

        SomUtils::MatD xEig(fullSz_, 1);
        RoptToEig(x, xEig);

        SomUtils::VecMatD R(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
        getRotations(xEig, R);

        SomUtils::MatD T(SomUtils::MatD::Zero(sz_.p_, sz_.n_));
        getTranslations(xEig, T);

        SomUtils::MatD Lambdas(SomUtils::MatD::Zero(numEdges_, 1));
        getScales(xEig, Lambdas);

        // SomUtils::MatD TijsScaled(SomUtils::MatD::Zero(sz_.d_, numEdges_));
        // ROFL_VAR1("makeTijsScaled")
        // makeTijsScaled(Tijs_, Lambdas, TijsScaled);
        // ROFL_VAR1(Tijs_);
        // ROFL_VAR1(TijsScaled);

        SomUtils::VecMatD rgR(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
        // ROFL_VAR3(xEig.rows(), R.size(), T.size());
        // ROFL_VAR3(R.size(), P.size(), rgR.size());
        rgradR(R, T, Lambdas, rgR);
        // ROFL_VAR5(rgR[0],rgR[1],rgR[2],rgR[3],rgR[4]);

        SomUtils::MatD rgT(SomUtils::MatD::Zero(sz_.p_, sz_.n_));
        rgradT(R, T, Lambdas, rgT);
        // ROFL_VAR1(rgT);

        SomUtils::MatD rgLambdas(SomUtils::MatD::Zero(numEdges_, 1));

        if (reluScaleCompensation_)
            rgradLambdasRelu(R, T, Lambdas, rgLambdas);
        else
            rgradLambdas(R, T, Lambdas, rgLambdas);

        // result->NewMemoryOnWrite();
        // result = Domain->RandomInManifold();

        int rotSz = getRotSz();
        int translSz = getTranslSz();

        int gElemIdx = 0;
        // fill result with computed gradient values: R
        for (int i = 0; i < sz_.n_; ++i)
        {
            // ROFL_VAR1(gElemIdx);
            // ROFL_VAR2("\n", rgR[gElemIdx]);
            // result->GetElement(gElemIdx).SetToIdentity(); // Ri
            // result->GetElement(gElemIdx).Print("Ri before assignment");

            Vector rgRiVec(sz_.p_, sz_.d_);
            // rgRiVec.Initialize();
            realdp *GroptlibWriteArray = rgRiVec.ObtainWriteEntireData();
            for (int j = 0; j < rotSz; ++j)
            {
                // ROFL_VAR2(i, j);
                // rgRiVec.Print("rgRiVec before assignment");

                // ROFL_VAR1(rgRiVec.GetElement(j, 0));

                GroptlibWriteArray[j] = rgR[i].reshaped(sz_.d_ * sz_.p_, 1)(j);

                // ROFL_VAR1("");
                // rgRiVec.Print("rgRiVec after assignment");
            }
            rgRiVec.CopyTo(result->GetElement(gElemIdx));
            // result->GetElement(gElemIdx).Print("Riem. grad Ri after assignment");
            gElemIdx++;
        }

        // fill result with computed gradient values: T

        Vector rgTiVec(sz_.p_, sz_.n_);
        realdp *GroptlibWriteArray = rgTiVec.ObtainWriteEntireData();
        for (int j = 0; j < sz_.p_ * sz_.n_; ++j)
        {
            // rgTiVec.Print("rgTiVec before assignment");

            // ROFL_VAR1(rgRiVec.GetElement(j, 0));

            GroptlibWriteArray[j] = rgT.reshaped(sz_.n_ * sz_.p_, 1)(j);

            // ROFL_VAR1("");
            // rgTiVec.Print("rgTiVec after assignment");
        }
        rgTiVec.CopyTo(result->GetElement(gElemIdx));
        gElemIdx++;
        // result->GetElement(gElemIdx).Print("grad Ti after assignment");

        // ROFL_VAR2("\n", rgT);

        // fill result with computed gradient values: Lambdas

        Vector rgLambdasIvec(numEdges_, 1);
        realdp *GroptlibWriteArray2 = rgLambdasIvec.ObtainWriteEntireData();
        for (int j = 0; j < numEdges_; ++j)
        {
            // rhTiVec.Print("rhTiVec before assignment");

            // ROFL_VAR1(rhRiVec.GetElement(j, 0));

            GroptlibWriteArray2[j] = rgLambdas(j); // TODO: reshaped() call can probably be removed

            // ROFL_VAR1("");
            // rhTiVec.Print("rhTiVec after assignment");
        }
        rgLambdasIvec.CopyTo(result->GetElement(gElemIdx));
        // result->GetElement(gElemIdx).Print("RieHess Lambdas after assignment");

        // ROFL_VAR1(gElemIdx);
        // result->GetElement(gElemIdx).Print();

        // result->NewMemoryOnWrite();
        // result->SetToZeros();
        // *result = Groptlib;

        // result->Print("printing final result");

        // ROFL_ASSERT(0);

        return *result;
    };

    // Vector &LsomProblem::Grad(const Variable &x, Vector *result) const
    // {
    // };

    void LsomProblem::makeTijsScaled(const SomUtils::MatD &tijs, const SomUtils::MatD &Lambdas, SomUtils::MatD &tijsScaled) const
    {
        ROFL_ASSERT_VAR5(tijs.rows() == tijsScaled.rows() && tijs.cols() == tijsScaled.cols() && Lambdas.rows() == tijsScaled.cols(), tijs.rows(), tijsScaled.rows(), tijs.cols(), tijsScaled.cols(), Lambdas.rows());
        tijsScaled = tijs;
        for (int e = 0; e < numEdges_; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1

            auto tij = tijs.col(e);
            double lambdaIJ = Lambdas(e, 0);

            tijsScaled.col(e) = tij * lambdaIJ;
        }
    }

    void LsomProblem::egradR(const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas,
                             SomUtils::VecMatD &egR) const
    {
        // num_edges = size(problem_data.edges, 1);
        // for
        //  e = 1 : num_edges
        //          ii = problem_data.edges(e, 1);
        //          jj = problem_data.edges(e, 2);
        //          Tj = T(:, jj);
        //          Ti = T(:, ii);
        //          lambdaij = lambdas(e, :);
        //          tij = problem_data.tijs( :, e);
        //          % R_i = R(:, :, ii);
        //          P_e = 2 * (Ti * lambdaij * tij ' - Tj * lambdaij * tij');
        //          g(:, :, ii) = g(:, :, ii) + P_e;
        for (int e = 0; e < numEdges_; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1

            auto Tj = T.col(j);
            auto Ti = T.col(i);
            double lambdaIJ = Lambdas(e, 0);

            auto tij = tijs_.col(e);

            auto P_e = 2 * lambdaIJ * (Ti * tij.transpose() - Tj * tij.transpose());

            egR[i] += P_e;
        }
    }

    void LsomProblem::rgradR(const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas,
                             SomUtils::VecMatD &rgR) const
    {
        SomUtils::VecMatD egR(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));

        egradR(R, T, Lambdas, egR);

        SomUtils::stiefelTangentProj(R, egR, rgR);
    }

    void LsomProblem::egradT(const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas,
                             SomUtils::MatD &egT) const
    {
        rgradT(R, T, Lambdas, egT);
    }

    void LsomProblem::rgradT(const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas,
                             SomUtils::MatD &rgT) const
    {
        // N = size(T, 2);
        // nrs = size(T, 1);

        // LR = zeros(N,N);
        // PR = zeros(nrs,N);
        // % BR_const = zeros(d,d);

        // num_edges = size(problem_data.edges,1);
        // for e = 1:num_edges
        //     ii = problem_data.edges(e,1);
        //     jj = problem_data.edges(e,2);
        //     bij = zeros(N,1);
        //     bij(ii, 1) = 1;
        //     bij(jj, 1) = -1;
        //     tij = problem_data.tijs(:, e);
        //     lambda_e = lambdas(e, 1);
        //     LR = LR + (bij * bij');
        //     Ri = R(:,:,ii);
        //     PR = PR + 2 * lambda_e * (Ri * tij * bij');

        // g=T*(LR+LR')+(PR);

        SomUtils::MatD LR = SomUtils::MatD::Zero(sz_.n_, sz_.n_);
        SomUtils::MatD PR = SomUtils::MatD::Zero(sz_.p_, sz_.n_);

        for (int e = 0; e < numEdges_; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1

            SomUtils::VecD bij = SomUtils::VecD::Zero(sz_.n_);
            bij(i) = 1.0;
            bij(j) = -1.0;

            auto tij = tijs_.col(e);
            double lambdaE = Lambdas(e, 0);
            LR += bij * bij.transpose();
            auto Ri = R[i];
            PR += 2 * lambdaE * (Ri * tij * bij.transpose());
        }
        rgT = T * (LR + LR.transpose()) + PR;
    }

    void LsomProblem::egradLambdasRelu(const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas,
                                       SomUtils::MatD &egLambdas) const
    {
        rgradLambdasRelu(R, T, Lambdas, egLambdas);
    }

    void LsomProblem::rgradLambdasRelu(const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas,
                                       SomUtils::MatD &rgLambdas) const
    {
        for (int e = 0; e < numEdges_; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1

            auto tij = tijs_.col(e);
            double lambdaIJ = Lambdas(e, 0);

            auto Ti = T.col(i);
            auto Tj = T.col(j);
            auto Ri = R[i];

            // base_part = 2*(tij_e'*tij_e * lambda_e + tij_e' * R_i' * T_i - tij_e' * R_i' * T_j);
            auto basePart = 2 * (tij.transpose() * tij * lambdaIJ + tij.transpose() * Ri.transpose() * Ti - tij.transpose() * Ri.transpose() * Tj);
            // ROFL_VAR1(basePart)

            // if lsom_relu_argument(lambda_e) > 0
            //     compensation_part = 2 * (lambda_e - 1);
            // else
            //     compensation_part = 0.0;
            double compensationPart = 0.0;
            if (lsomReLUargument(lambdaIJ) > 0)
            {
                // ROFL_VAR1("relu part active");
                // ROFL_VAR1(lambdaIJ);
                // ROFL_VAR1(lsomReLUargument(lambdaIJ));
                compensationPart = 2 * (lambdaIJ - 1.0);
            }
            // else
            // {
            // ROFL_VAR1("relu part inactive");
            // ROFL_VAR1(lambdaIJ);
            // ROFL_VAR1(lsomReLUargument(lambdaIJ));
            // double compensationPart = 0.0;
            // }

            // g_lambda(ee) = base_part + rho * compensation_part;
            rgLambdas(e, 0) = basePart(0, 0) + rho_ * compensationPart;
        }
        SomUtils::MatD lagrComp = -yAdmm_ + muAdmm_ * (Lambdas - zAdmm_); // Lagrange compensation
        rgLambdas += lagrComp;
    }

    void LsomProblem::egradLambdas(const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas,
                                   SomUtils::MatD &egLambdas) const
    {
        rgradLambdas(R, T, Lambdas, egLambdas);
    }

    void LsomProblem::rgradLambdas(const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas,
                                   SomUtils::MatD &rgLambdas) const
    {
        for (int e = 0; e < numEdges_; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1

            auto tij = tijs_.col(e);
            double lambdaE = Lambdas(e, 0);

            auto Ti = T.col(i);
            auto Tj = T.col(j);
            auto Ri = R[i];

            // base_part = 2*(tij_e'*tij_e * lambda_e + tij_e' * R_i' * T_i - tij_e' * R_i' * T_j);
            auto basePart = 2 * (tij.transpose() * tij * lambdaE + tij.transpose() * Ri.transpose() * Ti - tij.transpose() * Ri.transpose() * Tj);
            // ROFL_VAR1(basePart)
            double scaleCompensation = 0.0;

            // if l<1
            //     scale_compensation_ee= f=-1/(a*l-1)...
            // +1/(a-1)...
            // +b*(l-1);
            // else
            //     scale_compensation_ee=0;

            if (lambdaE <= 1 / a_)
            {
                if (fabs(rho_) < 1e-6)
                {
                    scaleCompensation = 0;
                }
                else
                {
                    scaleCompensation = std::nan("nan");
                }
            }
            else if (lambdaE <= 1.0)
            {
                scaleCompensation = (-1.0 / (a_ * lambdaE - 1.0)) + (1.0 / (a_ - 1.0)) + b_ * (lambdaE - 1.0);
            }

            // g_lambda(ee) = base_part + rho * compensation_part;
            rgLambdas(e, 0) = basePart(0, 0) + rho_ * scaleCompensation;
        }
        // lagrange_compensation = -y + mu * (lambdas - z);
        // g_lambda = g_lambda + lagrange_compensation;

        SomUtils::MatD lagrComp = -yAdmm_ + muAdmm_ * (Lambdas - zAdmm_); // Lagrange compensation
        rgLambdas += lagrComp;
    }

    double LsomProblem::lsomReLUargument(double lambdaE) const
    {
        return -lambdaE + 1.0;
    }

    void LsomProblem::computeHrt(const SomUtils::MatD &lambdas, const SomUtils::MatD &Tdot,
                                 SomUtils::VecMatD &h) const
    {
        // Ph = zeros(nrs, d * N);
        SomUtils::MatD Ph(SomUtils::MatD::Zero(sz_.p_, sz_.d_ * sz_.n_));

        // tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);
        // SomUtils::MatD TijsScaled(SomUtils::MatD::Zero(sz_.d_, numEdges_));
        // ROFL_VAR1("makeTijsScaled")
        // makeTijsScaled(Tijs_, lambdas, TijsScaled);

        // idx_col_p = reshape(1 : d * N, [], N)';

        // num_edges = size(problem_data.edges, 1);
        // for e = 1 : num_edges
        //     ii = problem_data.edges(e, 1);
        //     jj = problem_data.edges(e, 2);
        //     Tj_dot = Tdot(:, jj);
        //     Ti_dot = Tdot(:, ii);
        //     tij = tijs_scaled(:, e);
        //     P_e = 2 * (Ti_dot * tij' - Tj_dot * tij');
        //     Ph( :, idx_col_p(ii, :)) = Ph( :, idx_col_p(ii, :)) + P_e;

        int numEdges = edges_.rows();
        for (int e = 0; e < numEdges; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1

            auto Tj_dot = Tdot.col(j);
            auto Ti_dot = Tdot.col(i);
            auto tij = tijs_.col(e);
            double lambdaE = lambdas(e, 0);

            auto P_e = 2 * lambdaE * (Ti_dot - Tj_dot) * tij.transpose();
            // ROFL_VAR5(e, sz_.d_, sz_.p_, Ph.rows(), Ph.cols());
            // ROFL_VAR2(P_e.rows(), P_e.cols());
            Ph.block(0, i * sz_.d_, Ph.rows(), sz_.d_) += P_e;
        }
        SomUtils::unStackH(Ph, h, sz_.d_);
    }

    void LsomProblem::computeHrlambdas(const SomUtils::MatD &lambdasDot, const SomUtils::MatD &T,
                                       SomUtils::VecMatD &h) const
    {
        // Ph = zeros(nrs, d*N);
        SomUtils::MatD Ph(SomUtils::MatD::Zero(sz_.p_, sz_.d_ * sz_.n_));

        // tijs_dot_scaled = make_tijs_scaled(lambdas_dot, problem_data.tijs);
        // SomUtils::MatD TijsDotScaled(SomUtils::MatD::Zero(sz_.d_, numEdges_));
        // ROFL_VAR1("makeTijsScaled")
        // makeTijsScaled(Tijs_, lambdasDot, TijsDotScaled);

        // idx_col_p = reshape(1:d*N, [], N)';

        // num_edges = size(problem_data.edges,1);
        // for e = 1:num_edges
        //     ii = problem_data.edges(e,1);
        //     jj = problem_data.edges(e,2);
        //     T_j = T(:, jj);
        //     T_i = T(:, ii);
        //     tij_dot = tijs_dot_scaled(:,e);
        //     P_e = 2 * (T_i * tij_dot' - T_j * tij_dot');
        //     Ph(:, idx_col_p(ii, :)) = ...
        //         Ph(:, idx_col_p(ii, :)) + P_e;
        // end
        int numEdges = edges_.rows();
        for (int e = 0; e < numEdges; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1

            auto Tj = T.col(j);
            auto Ti = T.col(i);
            auto tij = tijs_.col(e);
            auto lambdaEdot = lambdasDot(e, 0);

            auto P_e = 2 * lambdaEdot * (Ti - Tj) * tij.transpose();
            Ph.block(0, i * sz_.d_, Ph.rows(), sz_.d_) += P_e;
        }

        SomUtils::unStackH(Ph, h, sz_.d_);
        // h=matUnstackH(Ph,d);
    }

    void LsomProblem::computeHtr(const SomUtils::MatD &lambdas, const SomUtils::VecMatD &uR,
                                 SomUtils::MatD &h) const
    {
        // tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);
        // [~, PR_dot] = make_LR_PR_BR_noloops(dR, tijs_scaled, problem_data.edges);
        // h=PR_dot';

        // SomUtils::MatD TijsScaled(SomUtils::MatD::Zero(sz_.d_, numEdges_));
        // ROFL_VAR1("makeTijsScaled")
        // makeTijsScaled(Tijs_, lambdas, TijsScaled);

        SomUtils::MatD PR_dot(SomUtils::MatD::Zero(sz_.n_, sz_.p_));

        for (int e = 0; e < numEdges_; ++e)
        {
            int ii = edges_(e, 0) - 1;
            int jj = edges_(e, 1) - 1;
            SomUtils::MatD uRi = uR[ii];

            SomUtils::MatD bij(SomUtils::MatD::Zero(sz_.n_, 1));

            bij(ii) = 1;
            bij(jj) = -1;

            double lambdaE = lambdas(e, 0);

            SomUtils::MatD tij = tijs_.col(e);

            PR_dot += 2 * lambdaE * bij * tij.transpose() * uRi.transpose(); // Matlab's tij is scaled!
        }

        h = PR_dot.transpose();
    }

    void LsomProblem::computeHtt(const SomUtils::MatD &lambdas, const SomUtils::VecMatD &xR, const SomUtils::MatD &uT,
                                 SomUtils::MatD &h) const
    {
        // tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);
        // [LR] = make_LR_PR_BR_noloops(R, tijs_scaled, problem_data.edges);
        // h = Tdot*(LR' + LR);

        // SomUtils::MatD TijsScaled(SomUtils::MatD::Zero(sz_.d_, numEdges_));
        // ROFL_VAR1("makeTijsScaled")
        // makeTijsScaled(Tijs_, lambdas, TijsScaled);

        // LR = zeros(N,N);
        // PR = zeros(N,nrs);
        // BR_const = zeros(d,d);
        SomUtils::MatD LR(SomUtils::MatD::Zero(sz_.n_, sz_.n_));

        for (int e = 0; e < numEdges_; ++e)
        {
            int ii = edges_(e, 0) - 1;
            int jj = edges_(e, 1) - 1;
            // SomUtils::MatD xRi = xR[ii];

            SomUtils::MatD bij(SomUtils::MatD::Zero(sz_.n_, 1));

            bij(ii) = 1;
            bij(jj) = -1;

            // double lambdaE = lambdas(e, 0);

            // SomUtils::MatD tij = tijs_.col(e);

            LR += bij * bij.transpose();
        }

        h = uT * (LR.transpose() + LR);
    }

    void LsomProblem::computeHtlambdas(const SomUtils::VecMatD &xR, const SomUtils::MatD &uLambdas,
                                       SomUtils::MatD &h) const
    {
        // h = zeros(size(T));
        // N = size(R,3);
        // num_edges = size(edges, 1);
        // for e = 1:num_edges
        //     ii = edges(e,1);
        //     jj = edges(e,2);
        //     Ri = R(:,:,ii);
        //     %         Rj = X.R(:,:,jj);
        //     %
        //     BIJ = zeros(N,1);
        //     BIJ(ii) = 1;
        //     BIJ(jj) = -1;
        //     %
        //     tij = problem_data.tijs(:, e);
        //     w_ij = BIJ * lambdas_dot(e) * tij' * Ri';
        //     h = h + 2 * w_ij';

        for (int e = 0; e < numEdges_; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1

            auto Ri = xR[i];

            SomUtils::MatD BIJ(SomUtils::MatD::Zero(sz_.n_, 1));
            BIJ(i) = 1;
            BIJ(j) = -1;

            auto tij = tijs_.col(e);
            auto w_ij = BIJ * uLambdas(e, 0) * tij.transpose() * Ri.transpose();
            h += 2 * w_ij.transpose();
        }
    }

    void LsomProblem::computeHlambdasr(const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR,
                                       const SomUtils::MatD &xT, const SomUtils::MatD &xLambdas,
                                       SomUtils::MatD &h) const
    {

        for (int e = 0; e < numEdges_; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1

            auto tij = tijs_.col(e);
            double lambdaE = xLambdas(e, 0);

            auto Ti = xT.col(i);
            auto Tj = xT.col(j);
            // auto Ri = xR[i];
            auto Ridot = uR[i];

            // e_th_elem_half = (tij' * R_i_dot') * (T_i - T_j) ;
            // eh(ee) = 2 * e_th_elem_half;

            double e_th_elem_half = (tij.transpose() * Ridot.transpose()) * (Ti - Tj); // 1x1 matrix

            h(e, 0) = 2 * e_th_elem_half;
        }
    }

    void LsomProblem::computeHlambdast(const SomUtils::VecMatD &xR, const SomUtils::MatD &uT,
                                       SomUtils::MatD &h) const
    {
        // for ee = 1:num_edges
        //     ii = edges(ee, 1);
        //     jj = edges(ee, 2);
        //     % lambdaE = x(ee);
        //     tij = tijs_vec(:, ee);
        //     % T_i = T(:,ii);
        //     % T_j = T(:,jj);
        //     T_i_dot = Tdot(:, ii);
        //     T_j_dot = Tdot(:, jj);
        //     % a = T_i - T_j;
        //     R_i = R(:, :, ii);
        //     b = R_i * tij;
        //     adot = T_i_dot - T_j_dot;
        //     e_th_elem = 2 * adot' * b;
        //     h(ee) = e_th_elem;

        for (int e = 0; e < numEdges_; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1

            auto tij = tijs_.col(e);
            auto Ti_dot = uT.col(i);
            auto Tj_dot = uT.col(j);

            auto Ri = xR[i];

            auto adot = Ti_dot - Tj_dot;
            auto b = Ri * tij;

            h(e, 0) = 2 * adot.transpose() * b;
        }
    }

    void LsomProblem::computeHlambdaslambdasRelu(const SomUtils::MatD &xLambdas, const SomUtils::MatD &uLambdas,
                                                 SomUtils::MatD &h) const
    {
        for (int e = 0; e < numEdges_; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            // int j = edges_(e, 1) - 1; // !! -1

            auto tij = tijs_.col(e);

            double lambdaE = xLambdas(e, 0);
            double lambdaDotE = uLambdas(e, 0);

            double compensationPart = 0.0;
            // if lsom_relu_argument(lambdas(ee)) > 0
            //     compensation_part = 2 * lambda_dot_ee;
            // else
            //     compensation_part = 0.0;
            if (lsomReLUargument(xLambdas(e, 0)) > 0)
            {
                // ROFL_VAR1("relu part active");
                // ROFL_VAR1(lambdaIJ);
                // ROFL_VAR1(lsomReLUargument(lambdaIJ));
                compensationPart = 2 * lambdaDotE;
            }

            // h(ee) = 2*lambda_dot_ee*(tij_e' * tij_e) + problem_data.rho * compensation_part;
            double basePart = 2 * lambdaDotE * (tij.transpose() * tij)(0, 0); // 1x1 matrix
            h(e, 0) = basePart + rho_ * compensationPart;
        }
        h += muAdmm_ * uLambdas;
    }

    void LsomProblem::computeHlambdaslambdas(const SomUtils::MatD &xLambdas, const SomUtils::MatD &uLambdas,
                                             SomUtils::MatD &h) const
    {
        for (int e = 0; e < numEdges_; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            // int j = edges_(e, 1) - 1; // !! -1

            auto tij = tijs_.col(e);

            double lambdaE = xLambdas(e, 0);
            double lambdaDotE = uLambdas(e, 0);

            double scaleCompensation = 0.0;
            // l = lambda_ee;
            // if l<=1
            //     compensation_part=a/(a*l-1)^2 + 0 + b;
            // else
            //     compensation_part=0;

            if (lambdaE <= 1 / a_)
            {
                if (fabs(rho_) < 1e-6)
                {
                    scaleCompensation = 0;
                }
                else
                {
                    scaleCompensation = std::nan("nan");
                }
            }
            else if (lambdaE <= 1.0)
            {
                scaleCompensation = a_ / ((a_ * lambdaE - 1.0) * (a_ * lambdaE - 1.0)) + b_;
            }

            // h(ee) = 2*lambda_dot_ee*(tij_e' * tij_e) + lambda_dot_ee * problem_data.rho * compensation_part;

            double basePart = 2 * lambdaDotE * (tij.transpose() * tij)(0, 0); // 1x1 matrix
            h(e, 0) = basePart + rho_ * lambdaDotE * scaleCompensation;
        }
        h += muAdmm_ * uLambdas;
    }

    void LsomProblem::hessGenprocEigen(const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR,
                                       const SomUtils::MatD &xT, const SomUtils::MatD &uT,
                                       const SomUtils::MatD &xLambdas, const SomUtils::MatD &uLambdas,
                                       SomUtils::VecMatD &rhR, SomUtils::MatD &rhT, SomUtils::MatD &rhLambdas) const
    {
        int staircaseStep = xT.rows();

        /*hrx*/
        // SomUtils::MatD hRR(SomUtils::MatD::Zero(staircaseStep, sz_.d_ * sz_.n_));
        // hRR = zeros(size(R));

        SomUtils::VecMatD hRT(sz_.n_, SomUtils::MatD::Zero(staircaseStep, sz_.d_));
        // ROFL_VAR1("computeHrt")
        computeHrt(xLambdas, uT, hRT);
        for (auto &hRTi : hRT)
        {
            // ROFL_VAR1(hRTi);
        }

        SomUtils::VecMatD hRLambdas(sz_.n_, SomUtils::MatD::Zero(staircaseStep, sz_.d_));
        // ROFL_VAR1("computeHrlambdas")
        computeHrlambdas(uLambdas, xT, hRLambdas);
        for (auto &hRlambdaI : hRLambdas)
        {
            // ROFL_VAR1(hRlambdaI);
        }

        /*htx*/
        SomUtils::MatD hTR(SomUtils::MatD::Zero(staircaseStep, sz_.n_));
        // ROFL_VAR1("computeHtr")
        computeHtr(xLambdas, uR, hTR);

        SomUtils::MatD hTT(SomUtils::MatD::Zero(staircaseStep, sz_.n_));
        // ROFL_VAR1("computeHtt")
        computeHtt(xLambdas, xR, uT, hTT);

        SomUtils::MatD hTLambdas(SomUtils::MatD::Zero(staircaseStep, sz_.n_));
        // ROFL_VAR1("computeHtlambdas")
        computeHtlambdas(xR, uLambdas, hTLambdas);

        /*hlambdasx*/
        SomUtils::MatD hLambdasR(SomUtils::MatD::Zero(numEdges_, 1));
        // ROFL_VAR1("computeHlambdasr")
        computeHlambdasr(xR, uR, xT, xLambdas, hLambdasR);

        SomUtils::MatD hLambdasT(SomUtils::MatD::Zero(numEdges_, 1));
        // ROFL_VAR1("computeHlambdasT")
        computeHlambdast(xR, uT, hLambdasT);

        SomUtils::MatD hLambdasLambdas(SomUtils::MatD::Zero(numEdges_, 1));
        // ROFL_VAR1("computeHlambdaslambdas")
        if (reluScaleCompensation_)
            computeHlambdaslambdasRelu(xLambdas, uLambdas, hLambdasLambdas);
        else
            computeHlambdaslambdas(xLambdas, uLambdas, hLambdasLambdas);

        // PUT EVERYTHING TOGETHER

        // ehR = lsom_ehess_R_R(R, Rdot, problem_data) + hrt + h_r_lambda;
        // egR = lsom_egrad_R(R, T, lambdas, problem_data);
        // h.R = manopt_stiefel_ehess2rhess(R, egR, ehR, Rdot);
        // h.T = lsom_ehess_T_T(R, T, Tdot, lambdas, problem_data) + htr + h_t_lambda;
        // h.lambda = h_lambda_lambda + h_lambda_r + h_lambda_t;
        auto ehR = SomUtils::VecMatD(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
        for (int i = 0; i < sz_.n_; ++i)
        {
            ehR[i] = hRT[i] + hRLambdas[i];
        }
        auto egR = SomUtils::VecMatD(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
        // ROFL_VAR1("lsomEgradR")
        lsomEgradR(xR, xT, xLambdas, egR);
        for (auto &egRi : egR)
        {
            // ROFL_VAR1(egRi);
        }
        for (auto &ehRi : ehR)
        {
            // ROFL_VAR1(ehRi);
        }

        manoptStiefelEhess2rhess(xR, egR, ehR, uR, rhR);
        for (auto &rhRi : rhR)
        {
            // ROFL_VAR1(rhRi);
        }

        rhT = hTR + hTT + hTLambdas;
        rhLambdas = hLambdasR + hLambdasT + hLambdasLambdas;
    }

    void LsomProblem::lsomEgradR(const SomUtils::VecMatD &xR, const SomUtils::MatD &xT,
                                 const SomUtils::MatD &xLambdas, SomUtils::VecMatD &egR) const
    {
        // P = zeros(nrs, d*N);
        SomUtils::MatD P(SomUtils::MatD::Zero(sz_.p_, sz_.d_ * sz_.n_));

        // tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);
        // SomUtils::MatD TijsScaled(SomUtils::MatD::Zero(sz_.d_, numEdges_));
        // ROFL_VAR1("lsom egradR")
        // ROFL_VAR1("makeTijsScaled")
        // makeTijsScaled(Tijs_, xLambdas, TijsScaled);

        // idx_col_p = reshape(1:d*N, [], N)';

        // num_edges = size(problem_data.edges,1);
        // for e = 1:num_edges
        //     ii = problem_data.edges(e,1);
        //     jj = problem_data.edges(e,2);
        //     T_j = T(:, jj);
        //     T_i = T(:, ii);
        //     tij = tijs_scaled(:,e);
        //     P_e = 2 * (T_i * tij' - T_j * tij');
        //     P(:, idx_col_p(ii, :)) = ...
        //         P(:, idx_col_p(ii, :)) + P_e;
        // end

        int numEdges = edges_.rows();
        for (int e = 0; e < numEdges; ++e)
        {
            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1

            auto Ti = xT.col(i);
            auto Tj = xT.col(j);
            auto tij = tijs_.col(e);
            auto lambdaE = xLambdas(e, 0);

            auto P_e = 2 * lambdaE * (Ti * tij.transpose() - Tj * tij.transpose());
            P.block(0, i * sz_.d_, P.rows(), sz_.d_) += P_e;
        }

        SomUtils::unStackH(P, egR, sz_.d_);
    }

    void LsomProblem::manoptStiefelEhess2rhess(const SomUtils::VecMatD &X,
                                               const SomUtils::VecMatD &egrad,
                                               const SomUtils::VecMatD &ehess,
                                               const SomUtils::VecMatD &Xdot,
                                               SomUtils::VecMatD &rH) const
    {
        // XtG = multiprod(multitransp(X), egrad);
        // symXtG = multisym(XtG);
        // HsymXtG = multiprod(H, symXtG);
        // rhess = stiefel_tangentProj(X, ehess - HsymXtG);
        for (int i = 0; i < sz_.n_; ++i)
        {
            manoptStiefelEhess2rhess(X[i], egrad[i], ehess[i], Xdot[i], rH[i]);
        }
    }

    void LsomProblem::manoptStiefelEhess2rhess(const SomUtils::MatD &X,
                                               const SomUtils::MatD &egrad,
                                               const SomUtils::MatD &ehess,
                                               const SomUtils::MatD &Xdot,
                                               SomUtils::MatD &rH) const
    {
        auto XtG = X.transpose() * egrad;
        auto symXtG = 0.5 * (XtG + XtG.transpose());
        auto HsymXtG = Xdot * symXtG;
        auto tmp = ehess - HsymXtG;
        SomUtils::stiefelTangentProj(X, tmp, rH);
    }

    void LsomProblem::hessGenprocEigenShifted(
        const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR,
        const SomUtils::MatD &xT, const SomUtils::MatD &uT,
        const SomUtils::MatD &xLambdas, const SomUtils::MatD &uLambdas,
        double mu,
        SomUtils::VecMatD &rhR, SomUtils::MatD &rhT, SomUtils::MatD &rhLambdas) const
    {
        hessGenprocEigen(xR, uR, xT, uT, xLambdas, uLambdas, rhR, rhT, rhLambdas);
        for (int i = 0; i < sz_.n_; ++i)
        {
            rhR[i] -= mu * uR[i];
        };
        // ROFL_VAR5(rhR[0], rhR[1], rhR[2], rhR[3], rhR[4]);

        // rhT = htt + htr;
        rhT -= mu * uT; // shift!
        // ROFL_VAR1(rhT);

        rhLambdas -= mu * uLambdas; // shift!
    }

    Vector &LsomProblem::RieHessianEta(const Variable &x, const Vector &etax, Vector *result) const
    {
        // TODO: implement

        result->NewMemoryOnWrite();

        // result->Print("RieHessianEta: printing it just after NewMemoryOnWrite()");

        // x.Print("x inside RieHessianEta");
        // etax.Print("etax inside RieHessianEta");

        SomUtils::MatD xEig(fullSz_, 1);
        RoptToEig(x, xEig);

        SomUtils::MatD xEtaEig(fullSz_, 1);
        RoptToEig(etax, xEtaEig);

        SomUtils::VecMatD R(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
        getRotations(xEig, R);

        SomUtils::VecMatD uR(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
        getRotations(xEtaEig, uR);

        SomUtils::MatD T(SomUtils::MatD::Zero(sz_.p_, sz_.n_));
        getTranslations(xEig, T);

        SomUtils::MatD uT(SomUtils::MatD::Zero(sz_.p_, sz_.n_));
        getTranslations(xEtaEig, uT);

        SomUtils::MatD Lambdas(SomUtils::MatD::Zero(numEdges_, 1));
        getScales(xEig, Lambdas);

        SomUtils::MatD uLambdas(SomUtils::MatD::Zero(numEdges_, 1));
        getScales(xEtaEig, uLambdas);

        // ROFL_VAR3(sz_.p_, xEig.rows(), xEig.cols())

        // ROFL_VAR2(R[0], uR[0]);
        // ROFL_VAR2(T, uT);

        SomUtils::VecMatD rhR(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
        SomUtils::MatD rhT(SomUtils::MatD::Zero(sz_.p_, sz_.n_));
        SomUtils::MatD rhLambdas(SomUtils::MatD::Zero(numEdges_, 1));
        hessGenprocEigen(R, uR, T, uT, Lambdas, uLambdas, rhR, rhT, rhLambdas);

        int rotSz = getRotSz();
        int translSz = getTranslSz();

        int gElemIdx = 0;
        // fill result with computed gradient values: R
        for (int i = 0; i < sz_.n_; ++i)
        {
            // ROFL_VAR1(gElemIdx);
            // ROFL_VAR2("\n", rhR[gElemIdx]);
            // result->GetElement(gElemIdx).SetToIdentity(); // Ri
            // result->GetElement(gElemIdx).Print("Ri before assignment");

            Vector rhRiVec(sz_.p_, sz_.d_);
            // rhRiVec.Initialize();
            realdp *GroptlibWriteArray = rhRiVec.ObtainWriteEntireData();
            for (int j = 0; j < rotSz; ++j)
            {
                // ROFL_VAR2(i, j);
                // rhRiVec.Print("rhRiVec before assignment");

                // ROFL_VAR1(rhRiVec.GetElement(j, 0));

                GroptlibWriteArray[j] = rhR[i].reshaped(sz_.d_ * sz_.p_, 1)(j);

                // ROFL_VAR1("");
                // rhRiVec.Print("rhRiVec after assignment");
            }
            rhRiVec.CopyTo(result->GetElement(gElemIdx));
            // result->GetElement(gElemIdx).Print("Riem. grad Ri after assignment");
            gElemIdx++;
        }

        // fill result with computed gradient values: T

        Vector rhTiVec(sz_.p_, sz_.n_);
        realdp *GroptlibWriteArray = rhTiVec.ObtainWriteEntireData();
        for (int j = 0; j < sz_.p_ * sz_.n_; ++j)
        {
            // rhTiVec.Print("rhTiVec before assignment");

            // ROFL_VAR1(rhRiVec.GetElement(j, 0));

            GroptlibWriteArray[j] = rhT.reshaped(sz_.n_ * sz_.p_, 1)(j);

            // ROFL_VAR1("");
            // rhTiVec.Print("rhTiVec after assignment");
        }
        rhTiVec.CopyTo(result->GetElement(gElemIdx));
        gElemIdx++;
        // result->GetElement(gElemIdx).Print("RieHess T after assignment");

        // fill result with computed gradient values: Lambdas

        Vector rhLambdasIvec(numEdges_, 1);
        realdp *GroptlibWriteArray2 = rhLambdasIvec.ObtainWriteEntireData();
        for (int j = 0; j < numEdges_; ++j)
        {
            // rhTiVec.Print("rhTiVec before assignment");

            // ROFL_VAR1(rhRiVec.GetElement(j, 0));

            GroptlibWriteArray2[j] = rhLambdas.reshaped(numEdges_, 1)(j); // TODO: reshaped() call can probably be removed

            // ROFL_VAR1("");
            // rhTiVec.Print("rhTiVec after assignment");
        }
        rhLambdasIvec.CopyTo(result->GetElement(gElemIdx));
        // result->GetElement(gElemIdx).Print("RieHess Lambdas after assignment");

        // result->Print("RieHessianEta: printing just before end of function");

        return *result;
    };

    void LsomProblem::RoptToEig(Vector x, SomUtils::MatD &xEigen) const
    {
        Vector xT = x.GetTranspose(); // Eigen ADV init is row-major!!

        int totSz = xEigen.rows();

        const realdp *xArr = xT.ObtainWriteEntireData();
        for (int i = 0; i < totSz; ++i)
            xEigen(i) = xArr[i];
    }

    void LsomProblem::getRi(const Variable &x, SomUtils::MatD &rOut, int i) const
    {
        SomUtils::MatD xEigen(fullSz_, 1);
        RoptToEig(x, xEigen);

        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // as a vector

        int startId = i * rotSz;
        // int endId = (i+1) * rotSz;

        SomUtils::MatD rOutVec(SomUtils::MatD::Zero(rotSz, 1));
        rOutVec = xEigen.block(startId, 0, rotSz, 1);

        rOut = rOutVec.reshaped(sz_.p_, sz_.d_);
    }

    void LsomProblem::getRi(const SomUtils::MatD &xEig, SomUtils::MatD &rOut, int i) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // as a vector

        int startId = i * rotSz;
        // int endId = (i+1) * rotSz;

        SomUtils::MatD rOutVec(SomUtils::MatD::Zero(rotSz, 1));
        rOutVec = xEig.block(startId, 0, rotSz, 1);

        rOut = rOutVec.reshaped(sz_.p_, sz_.d_);
    }

    void LsomProblem::getRiSEdN(const SomUtils::MatD &xEig, SomUtils::MatD &rOut, int i) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = sz_.d_ * sz_.d_; // as a vector

        int startId = i * rotSz;
        // int endId = (i+1) * rotSz;

        SomUtils::MatD rOutVec(SomUtils::MatD::Zero(rotSz, 1));
        rOutVec = xEig.block(startId, 0, rotSz, 1);

        rOut = rOutVec.reshaped(sz_.d_, sz_.d_);
    }

    void LsomProblem::getRotations(const SomUtils::MatD &xEig, SomUtils::VecMatD &rOut) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // vectorized

        ROFL_ASSERT_VAR5(rotSz * sz_.n_ + sz_.p_ * sz_.n_ + numEdges_ == xEig.rows(),
                         rotSz, sz_.n_, sz_.p_, numEdges_, xEig.rows());

        // int endId = (i+1) * rotSz;

        std::for_each(rOut.begin(), rOut.end(), [](SomUtils::MatD &x) { //^^^ take argument by reference: LAMBDA FUNCTION
            x.setZero();
        });

        for (int i = 0; i < sz_.n_; ++i)
        {
            int startId = i * rotSz;
            SomUtils::MatD Ri(SomUtils::MatD::Zero(sz_.p_, sz_.d_));
            getRi(xEig, Ri, i); // TODO: this can probably be optimized better
            rOut[i] = Ri;
            // ROFL_VAR3(i, rOut[i].rows(), rOut[i].cols())
        }
    }

    void LsomProblem::getTi(const Variable &x, SomUtils::MatD &tOut, int i) const
    {
        SomUtils::MatD xEigen(fullSz_, 1);
        RoptToEig(x, xEigen);

        // rOut already needs to have fixed size by here
        int rotSz = getRotSz();
        int translSz = getTranslSz();

        int startId = sz_.n_ * rotSz + i * translSz;
        // int endId = (i+1) * rotSz;

        // ROFL_VAR1(startId);

        tOut.setZero(); // TODO: remove later
        tOut = xEigen.block(startId, 0, translSz, 1);
    }

    void LsomProblem::getTi(const SomUtils::MatD &xEig, SomUtils::MatD &tOut, int i) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz();
        int translSz = getTranslSz();

        int startId = sz_.n_ * rotSz + i * translSz;
        // int endId = (i+1) * rotSz;

        // ROFL_VAR1(startId);
        tOut.setZero(); // TODO: remove later
        tOut = xEig.block(startId, 0, translSz, 1);
    }

    void LsomProblem::getLambdaI(const SomUtils::MatD &xEig, double &lambdaOut, int i) const
    {
        int rotSz = getRotSz();
        int translSz = getTranslSz();

        int id = sz_.n_ * rotSz + sz_.n_ * translSz + i;

        // ROFL_VAR2(id, xEig.rows())
        lambdaOut = xEig(id, 0);
        // ROFL_VAR1(lambdaOut);
    }

    void LsomProblem::getTiSEdN(const SomUtils::MatD &xEig, SomUtils::MatD &tOut, int i) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = sz_.d_ * sz_.d_;
        int translSz = sz_.d_;

        int startId = sz_.n_ * rotSz + i * translSz;
        // int endId = (i+1) * rotSz;

        // ROFL_VAR1(startId);
        tOut.setZero(); // TODO: remove later
        tOut = xEig.block(startId, 0, translSz, 1);
    }

    void LsomProblem::getLambdaISEdN(const SomUtils::MatD &xEig, double &lambdaOut, int i) const
    {
        int rotSz = sz_.d_ * sz_.d_;
        int translSz = sz_.d_;

        int id = sz_.n_ * rotSz + sz_.n_ * translSz + i + 1;

        lambdaOut = xEig(id, 0);
        // ROFL_VAR1(lambdaOut);
    }

    void LsomProblem::getTranslations(const SomUtils::MatD &xEig, SomUtils::MatD &tOut) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // as a vector

        // int endId = (i+1) * rotSz;

        for (int i = 0; i < sz_.n_; ++i)
        {
            SomUtils::MatD Ti(sz_.p_, 1);
            getTi(xEig, Ti, i); // TODO: this can probably be optimized better
            tOut.col(i) = Ti;
        }
    }

    void LsomProblem::getScales(const SomUtils::MatD &xEig, SomUtils::MatD &scalesOut) const
    {
        // rOut already needs to have fixed size by here
        // int rotSz = getRotSz();       // as a vector
        // int translSz = getTranslSz(); // as a vector

        // int endId = (i+1) * rotSz;

        for (int i = 0; i < numEdges_; ++i)
        {
            double scaleI = 0.0;
            getLambdaI(xEig, scaleI, i); // TODO: this can probably be optimized better
            scalesOut(i, 0) = scaleI;
        }
    }

    void LsomProblem::makePfrct(const SomUtils::MatD &T, const SomUtils::MatD &Lambdas,
                                SomUtils::MatD &P, double &frct) const
    {
        // P = zeros(nrs, d*N);

        // TODO: add some size checks at least for P

        P.setZero();
        frct = 0.0;

        for (int e = 0; e < numEdges_; ++e)
        {
            int ii = edges_(e, 0) - 1;
            int jj = edges_(e, 1) - 1;

            SomUtils::MatD T_j = T.col(jj);
            SomUtils::MatD T_i = T.col(ii);
            SomUtils::MatD tij = tijs_.col(e);
            double lambdaE = Lambdas(e, 0);
            SomUtils::MatD Pe = 2 * lambdaE * (T_i * tij.transpose() - T_j * tij.transpose());
            P.block(0, ii * sz_.d_, T.rows(), sz_.d_) += Pe;

            SomUtils::MatD c = T_i * T_i.transpose() + T_j * T_j.transpose() - T_i * T_j.transpose() - T_j * T_i.transpose();
            SomUtils::MatD d = lambdaE * lambdaE * tij * tij.transpose();
            frct += c.trace() + d.trace();
        }
    }

    void LsomProblem::makeLrPrBr(const SomUtils::VecMatD &R, const SomUtils::MatD &Lambdas,
                                 SomUtils::MatD &Lr, SomUtils::MatD &Pr, SomUtils::MatD &Br) const
    {
        // LR = zeros(N,N);
        // PR = zeros(N,nrs);
        // BR_const = zeros(d,d);

        // TODO: add some size checks at least for Lr, Pr, Br

        Lr.setZero();
        Pr.setZero();
        Br.setZero();

        for (int e = 0; e < numEdges_; ++e)
        {
            int ii = edges_(e, 0) - 1;
            int jj = edges_(e, 1) - 1;
            SomUtils::MatD Ri = R[ii];

            SomUtils::MatD bij(SomUtils::MatD::Zero(sz_.n_, 1));

            bij(ii) = 1;
            bij(jj) = -1;

            double lambdaE = Lambdas(e, 0);

            SomUtils::MatD tij = tijs_.col(e);

            Lr += bij * bij.transpose();

            Pr += 2 * lambdaE * Ri * tij * bij.transpose();

            Br += lambdaE * lambdaE * tij * tij.transpose();
        }
    }

    int LsomProblem::getRotSz() const
    {
        return sz_.d_ * sz_.p_;
    }

    int LsomProblem::getTranslSz() const
    {
        return sz_.p_;
    }

    void LsomProblem::setRho(double rho)
    {
        rho_ = rho;
    }

    void LsomProblem::setUsePIM(bool usePIM)
    {
        usePIM_ = usePIM;
    }

    void LsomProblem::setPimMaxIterations(int numMaxIterations)
    {
        pimMaxIterations_ = numMaxIterations;
    }

    void LsomProblem::setReluScaleCompensation(bool flag)
    {
        reluScaleCompensation_ = flag;
    }

    void LsomProblem::setLogScaleCompensationParam(double a)
    {
        a_ = a;
        b_ = -a_ / ((a_ - 1) * (a_ - 1)); // dependent on a_
    }

    void LsomProblem::setZAdmm(const SomUtils::MatD &z)
    {
        zAdmm_ = z;
    }

    void LsomProblem::setYAdmm(const SomUtils::MatD &y)
    {
        yAdmm_ = y;
    }

    void LsomProblem::setMuAdmm(double mu)
    {
        muAdmm_ = mu;
    }

    void LsomProblem::setMaxIterAdmm(int maxIterAdmm)
    {
        maxIterAdmm_ = maxIterAdmm;
    }

    void LsomProblem::setPerformGlobalization(bool performGlobalization)
    {
        performGlobalization_ = performGlobalization;
    }

    void LsomProblem::setSsomInitguess(bool ssomInitguess)
    {
        ssomInitguess_ = ssomInitguess;
    }

    void LsomProblem::setFirstZadmmLambdas(bool firstZadmmLambdas)
    {
        firstZadmmLambdas_ = firstZadmmLambdas;
    }

    void LsomProblem::vectorizeR(const SomUtils::VecMatD &R, SomUtils::MatD &RvecOut) const
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

    void LsomProblem::vectorizeRTLambdas(const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas,
                                         SomUtils::MatD &XvecOut) const
    {
        // int fullRotsSz = sz_.p_ * sz_.d_ * sz_.n_;

        // for (int i=0; i<fullRotsSz; ++i) {
        // }

        int n = R.size();

        ROFL_ASSERT(T.cols() == n)

        int fullIdx = 0;
        for (int i = 0; i < n; ++i)
        {
            int ric = R[i].cols();
            int rir = R[i].rows();
            for (int j = 0; j < ric; ++j)
            {
                for (int k = 0; k < rir; ++k)
                {
                    XvecOut(fullIdx, 0) = R[i](k, j);
                    fullIdx++;
                    // ROFL_VAR4(i, j, k, fullIdx);
                }
            }
        }

        for (int i = 0; i < T.cols(); ++i)
        {
            for (int j = 0; j < T.rows(); ++j)
            {
                XvecOut(fullIdx, 0) = T(j, i);
                fullIdx++;
                // ROFL_VAR4(i, j, k, fullIdx);
            }
        }

        for (int i = 0; i < Lambdas.rows(); ++i)
        {
            XvecOut(fullIdx, 0) = Lambdas(i, 0);
            fullIdx++;
            // ROFL_VAR4(i, j, k, fullIdx);
        }

        // TODO: more asserts may be added

        ROFL_ASSERT(fullIdx == XvecOut.rows())
    }

    void LsomProblem::RoptToEigStiefel(Vector x, SomUtils::MatD &xEigen) const
    {
        Vector xT = x.GetTranspose(); // Eigen ADV init is row-major!!

        int totSz = xEigen.rows(); // xEigen is supposed to be a vectorized matrix

        const realdp *xArr = xT.ObtainWriteEntireData();
        for (int i = 0; i < totSz; ++i)
            xEigen(i) = xArr[i];
    }

    void LsomProblem::makeAdjMatFromEdges(Eigen::MatrixXi &adjMat) const
    {
        ROFL_ASSERT(adjMat.rows() == sz_.n_)
        ROFL_ASSERT(adjMat.cols() == sz_.n_)

        adjMat.setZero();
        for (int k = 0; k < numEdges_; ++k)
        {
            int ii = edges_(k, 0) - 1;
            int jj = edges_(k, 1) - 1;
            adjMat(ii, jj) = 1;
        }
    }

    void LsomProblem::setGtR(const SomUtils::VecMatD &R)
    {
        Rgt_ = R;
    }

    void LsomProblem::setGtT(const SomUtils::MatD &T)
    {
        Tgt_ = T;
    }

    void LsomProblem::setGtLambdas(const SomUtils::MatD &Lambdas)
    {
        LambdasGt_ = Lambdas;
    }

    void LsomProblem::setGt(const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas)
    {
        Rgt_ = R;
        Tgt_ = T;
        LambdasGt_ = Lambdas;
    }

    bool LsomProblem::getRsRecoverySuccess() const
    {
        return rsRecoverySuccess_;
    }

    void LsomProblem::setCostCurr(double cc)
    {
        costCurr_ = cc;
    }

    void LsomProblem::setEnableRs(bool enableRs)
    {
        enableRs_ = enableRs;
    }

    void LsomProblem::setTolerancesAdmm(double tolPrimal, double tolDual)
    {
        tolAdmmPrimal_ = tolPrimal;
        tolAdmmDual_ = tolDual;
    }

    void LsomProblem::makeHmat(const SomUtils::MatD &XvecNext, const SomUtils::SomSize &szNext, SomUtils::MatD &Hmat) const
    {
        int staircaseStepLevel = szNext.p_; // TODO: it can maybe be deducted fron XvecNext.size()?

        ROFL_VAR1(staircaseStepLevel)

        ROPTLIB::LsomProblem ProbNextLocal(szNext, tijs_, edges_);

        SomUtils::VecMatD xR(szNext.n_, SomUtils::MatD::Zero(staircaseStepLevel, szNext.d_));
        SomUtils::MatD xT(SomUtils::MatD::Zero(staircaseStepLevel, szNext.n_));
        SomUtils::MatD xLambdas(SomUtils::MatD::Zero(numEdges_, 1));
        ROFL_VAR1("Calling getRotations()")
        ProbNextLocal.getRotations(XvecNext, xR);
        ProbNextLocal.getTranslations(XvecNext, xT);
        ProbNextLocal.getScales(XvecNext, xLambdas);

        int vecsz = XvecNext.rows();
        ROFL_VAR1(vecsz)

        for (int i = 0; i < vecsz; ++i)
        {
            SomUtils::MatD eI(SomUtils::MatD::Zero(vecsz, 1));
            eI(i) = 1;
            SomUtils::VecMatD uRi(szNext.n_, SomUtils::MatD::Zero(staircaseStepLevel, szNext.d_));
            SomUtils::MatD uTi(SomUtils::MatD::Zero(staircaseStepLevel, szNext.n_));
            SomUtils::MatD uLambdas(SomUtils::MatD::Zero(numEdges_, 1));

            ROFL_VAR1("Calling getRotations()")
            ProbNextLocal.getRotations(eI, uRi);
            ProbNextLocal.getTranslations(eI, uTi);
            ProbNextLocal.getScales(eI, uLambdas);

            SomUtils::VecMatD rhrI(szNext.n_, SomUtils::MatD::Zero(staircaseStepLevel, szNext.d_));
            SomUtils::MatD rhtI(SomUtils::MatD::Zero(staircaseStepLevel, szNext.n_));
            SomUtils::MatD rhLambdasI(SomUtils::MatD::Zero(numEdges_, 1));

            ROFL_VAR2(i, vecsz)
            ROFL_VAR1(tijs_)
            ROFL_VAR3(szNext.n_, xR[szNext.n_ - 1].rows(), xR[szNext.n_ - 1].cols())
            ROFL_VAR2(uRi[0], uTi)
            ProbNextLocal.hessGenprocEigen(xR, uRi, xT, uTi, xLambdas, uLambdas, rhrI, rhtI, rhLambdasI);
            ROFL_VAR2(rhrI[0], rhtI)

            SomUtils::MatD rhVecI(SomUtils::MatD::Zero(vecsz, 1));
            ProbNextLocal.vectorizeRTLambdas(rhrI, rhtI, rhLambdasI, rhVecI);
            Hmat.col(i) = rhVecI;
            ROFL_VAR2(i, rhVecI.transpose())
        }
    }

    void LsomProblem::makeHmatLsom(const SomUtils::MatD &XvecNext, const SomUtils::SomSize &szNext, SomUtils::MatD &Hmat) const
    {
        int staircaseStepLevel = szNext.p_; // TODO: it can maybe be deducted fron XvecNext.size()?

        ROFL_VAR1(staircaseStepLevel)

        ROPTLIB::LsomProblem ProbNextLocal(szNext, tijs_, edges_);

        SomUtils::VecMatD xR(szNext.n_, SomUtils::MatD::Zero(staircaseStepLevel, szNext.d_));
        SomUtils::MatD xT(SomUtils::MatD::Zero(staircaseStepLevel, szNext.n_));
        SomUtils::MatD xLambdas(SomUtils::MatD::Zero(numEdges_, 1));
        ROFL_VAR1("Calling getRotations()")
        ProbNextLocal.getRotations(XvecNext, xR);
        ProbNextLocal.getTranslations(XvecNext, xT);
        ProbNextLocal.getScales(XvecNext, xLambdas);

        int vecsz = XvecNext.rows();
        ROFL_VAR1(vecsz)

        // for (int i = 0; i < vecsz; ++i)
        // {
        //     SomUtils::MatD eI(SomUtils::MatD::Zero(vecsz, 1));
        //     eI(i) = 1;
        //     SomUtils::VecMatD uRi(szNext.n_, SomUtils::MatD::Zero(staircaseStepLevel, szNext.d_));
        //     SomUtils::MatD uTi(SomUtils::MatD::Zero(staircaseStepLevel, szNext.n_));
        //     SomUtils::MatD uLambdas(SomUtils::MatD::Zero(numEdges_, 1));

        //     ROFL_VAR1("Calling getRotations()")
        //     ProbNextLocal.getRotations(eI, uRi);
        //     ProbNextLocal.getTranslations(eI, uTi);
        //     ProbNextLocal.getScales(eI, uLambdas);

        //     SomUtils::VecMatD rhrI(szNext.n_, SomUtils::MatD::Zero(staircaseStepLevel, szNext.d_));
        //     SomUtils::MatD rhtI(SomUtils::MatD::Zero(staircaseStepLevel, szNext.n_));
        //     SomUtils::MatD rhLambdasI(SomUtils::MatD::Zero(numEdges_, 1));

        //     ROFL_VAR2(i, vecsz)
        //     ROFL_VAR1(tijs_)
        //     ROFL_VAR3(szNext.n_, xR[szNext.n_ - 1].rows(), xR[szNext.n_ - 1].cols())
        //     ROFL_VAR2(uRi[0], uTi)
        //     ProbNextLocal.hessGenprocEigen(xR, uRi, xT, uTi, xLambdas, uLambdas, rhrI, rhtI, rhLambdasI);
        //     ROFL_VAR2(rhrI[0], rhtI)

        //     SomUtils::MatD rhVecI(SomUtils::MatD::Zero(vecsz, 1));
        //     ProbNextLocal.vectorizeRTLambdas(rhrI, rhtI, rhLambdasI, rhVecI);
        //     Hmat.col(i) = rhVecI;
        //     ROFL_VAR2(i, rhVecI.transpose())
        // }
    }

    void LsomProblem::updateLsomPenaltyParam(const double &mu, const SomUtils::MatD &xK, const SomUtils::MatD &zK, const SomUtils::MatD &zPrev,
                                             double &muNext, SomUtils::MatD &rK, SomUtils::MatD &sK) const
    {
        // r_k = x_k - z_k;
        // s_k = -mu_prev * (z_k - z_prev);
        rK = xK - zK;
        sK = -mu * (zK - zPrev);

        // tau_incr = 2.0;
        // tau_decr = 2.0;
        // mu_tau = 1.0;

        double tauIncr = 2.0; // TODO: make these settable from outside
        double tauDecr = 2.0; // TODO: make these settable from outside
        double muTau = 1.0;   // TODO: make these settable from outside

        // assert(~ ((norm(r_k) > mu_tau * norm(s_k)) && (norm(s_k) > mu_tau * norm(r_k))) )
        ROFL_ASSERT_VAR3(!(rK.norm() > muTau * sK.norm() && sK.norm() > muTau * rK.norm()), rK.norm(), sK.norm(), muTau)

        // if norm(r_k) > mu_tau * norm(s_k)
        //     mu_next = tau_incr * mu_prev;
        // elseif norm(s_k) > mu_tau * norm(r_k)
        //     mu_next = mu_prev / tau_decr;
        // else
        //     mu_next = mu_prev;

        if (rK.norm() > muTau * sK.norm())
        {
            muNext = tauIncr * mu;
        }
        else if (sK.norm() > muTau * rK.norm())
        {
            muNext = mu / tauDecr;
        }
        else
        {
            muNext = mu;
        }

        // % mu_next = max(0.1, mu_next); % This puts a cap on mu minimum value

        // % mu_next = mu_prev; % TODO: this clears mu_next updates
    }

    bool LsomProblem::checkAdmmStoppingCondition(const SomUtils::MatD &rK, const SomUtils::MatD &sK, double epsAbs, double epsRel) const
    {
        // A = eye(num_edges);
        // B = eye(num_edges);
        // % C = zeros(num_edges, 1); % excluding norm(C) from further consideration
        SomUtils::MatD A = SomUtils::MatD::Identity(numEdges_, numEdges_);
        SomUtils::MatD B = SomUtils::MatD::Identity(numEdges_, numEdges_);

        // tmp = [norm(A * x_k), norm(B*z_k)];
        SomUtils::MatD tmp = SomUtils::MatD::Zero(1, 2);
        tmp(0, 0) = (A * rK).norm();
        tmp(0, 1) = (B * sK).norm();
        // eps_pri = sqrt(num_edges) * eps_abs + eps_rel * max(tmp, [], "all");
        // eps_dual = sqrt(num_edges) * eps_abs + eps_rel * norm(A' * y_k);
        double epsPri = std::sqrt(numEdges_) * epsAbs + epsRel * tmp.maxCoeff();
        double epsDual = std::sqrt(numEdges_) * epsAbs + epsRel * (A.transpose() * sK).norm();

        // stop = norm(r_k) <= eps_pri && norm(s_k) <= eps_dual;
        bool stop = rK.norm() <= epsPri && sK.norm() <= epsDual;

        return stop;
    }

    double runLsom(ROPTLIB::LsomProblem &Prob,
                   const ROPTLIB::Vector &startX,
                   int src,
                   SomUtils::VecMatD &Rout,
                   SomUtils::MatD &Tout,
                   SomUtils::MatD &lambdasOut,
                   int &staircaseStepIdx,
                   bool &rsSuccess, bool &rotDetsOk, bool &lambdasAcceptable, bool &rsActuallyUseful)
    {
        ROFL_VAR1("Start of runLsom()")

        if (Prob.reluScaleCompensation_)
        {
            ROFL_VAR1(Prob.costEigenRelu(Prob.Rgt_, Prob.Tgt_, Prob.LambdasGt_));
        }
        else
        {
            ROFL_VAR1(Prob.costEigen(Prob.Rgt_, Prob.Tgt_, Prob.LambdasGt_));
        }

        auto rGt = Prob.Rgt_;
        auto tGt = Prob.Tgt_;
        auto lambdasGt = Prob.LambdasGt_;

        int d = Prob.sz_.d_;
        int n = Prob.sz_.n_;
        int e = Prob.numEdges_;

        double costOut = std::numeric_limits<double>::infinity();

        bool admmStoppingConditionReached = false;

        SomUtils::MatD zAdmm(SomUtils::MatD::Zero(e, 1));

        // zAdmm = Prob.zAdmm_;
        if (Prob.firstZadmmLambdas_)
        {
            SomUtils::MatD startXEig(Prob.fullSz_, 1);
            Prob.RoptToEig(startX, startXEig);
            SomUtils::MatD lambdasInitguess = SomUtils::MatD::Zero(e, 1);
            Prob.getScales(startXEig, lambdasInitguess);
            zAdmm = lambdasInitguess; // TODO: use setter instead of direct access
        }
        else
        {
            zAdmm = Prob.zAdmm_; // TODO: use setter instead of direct access
        }

        SomUtils::MatD yAdmm(SomUtils::MatD::Zero(e, 1));
        yAdmm = Prob.yAdmm_;
        double muAdmm = Prob.muAdmm_;

        ROFL_VAR1("ADMM init:")
        ROFL_VAR1(zAdmm.transpose())
        ROFL_VAR1(yAdmm.transpose())
        ROFL_VAR1(muAdmm)

        ROPTLIB::Vector startXlocal = startX;
        // startX.CopyTo(startXlocal);

        if (Prob.ssomInitguess_)
        {
            SsomProblem ssomProb(Prob.sz_, Prob.tijs_, Prob.edges_);
            ssomProb.setReluScaleCompensation(Prob.reluScaleCompensation_);

            // integer numoftypes = 3;
            // ROPTLIB::Stiefel mani1next(d, d);
            // mani1next.ChooseParamsSet2();
            // ROPTLIB::Euclidean mani2next(somSzNext.p_, somSzNext.n_);
            // ROPTLIB::ProductManifold ProdMani(numoftypes,
            //                                       &mani1next, numofmani1, &mani2next, numofmani2, &mani3, numofmani3);
            ssomProb.SetDomain(Prob.GetDomain());

            ssomProb.setUsePIM(true);           // same as default
            ssomProb.setPimMaxIterations(5000); // same as default
            ssomProb.setGt(rGt, tGt, lambdasGt);

            ssomProb.setRho(0.0); // TODO: make this settable from outside, and maybe also make it adaptive as in ADMM?

            ROFL_VAR1("Using SSOM initguess")

            if (ssomProb.reluScaleCompensation_)
            {
                ROFL_VAR1(ssomProb.costEigenRelu(ssomProb.Rgt_, ssomProb.Tgt_, ssomProb.LambdasGt_));
            }
            else
            {
                ROFL_VAR1(ssomProb.costEigen(ssomProb.Rgt_, ssomProb.Tgt_, ssomProb.LambdasGt_));
            }

            // output the parameters of the manifold of domain
            startXlocal.Print("startXlocal before SSOM optimization");
            ROPTLIB::RTRNewton *RTRNewtonSolver = new ROPTLIB::RTRNewton(&ssomProb, &startXlocal);
            RTRNewtonSolver->Verbose = ROPTLIB::ITERRESULT;
            // RTRNewtonSolver->Max_Iteration = 500;
            // RTRNewtonSolver->Max_Inner_Iter = 500;
            // ROPTLIB::PARAMSMAP solverParams = {std::pair<std::string, double>("Max_Inner_Iter", 10)};
            // RTRNewtonSolver->SetParams(solverParams);
            RTRNewtonSolver->CheckParams();

            RTRNewtonSolver->Run();
            // Numerically check gradient consistency (optional).
            auto Xopt = RTRNewtonSolver->GetXopt();
            auto XoptCost = RTRNewtonSolver->Getfinalfun();

            Xopt.CopyTo(startXlocal);

            delete RTRNewtonSolver;
        }

        auto Xbackup = startXlocal;

        if (!Prob.enableRs_)
        {
            std::cout << "RS disabled: skipping staircase and returning directly with costOut = costLast" << std::endl;

            int iterAdmm = 0;

            SomUtils::MatD rK(SomUtils::MatD::Zero(e, 1));
            SomUtils::MatD sK(SomUtils::MatD::Zero(e, 1));

            // while iter_admm < 100 && ~admm_stopping_condition_reached
            while (iterAdmm < Prob.maxIterAdmm_ && !admmStoppingConditionReached)
            {
                ROFL_VAR2(iterAdmm, Prob.maxIterAdmm_)

                auto ProbAdmm = Prob;
                ProbAdmm.costCurr_ = costOut; // TODO: use setter instead of direct access
                ProbAdmm.setZAdmm(zAdmm);
                ProbAdmm.setYAdmm(yAdmm);
                ProbAdmm.setMuAdmm(muAdmm);

                // lsomRTR(nrs, d, N, problem_data, params, transf_initguess_struct, lambdas_initguess)
                ROPTLIB::RTRNewton *RTRNewtonSolver = new ROPTLIB::RTRNewton(&ProbAdmm, &startXlocal); // USE INITGUESS HERE!
                RTRNewtonSolver->Verbose = ROPTLIB::ITERRESULT;
                // RTRNewtonSolver->Max_Iteration = 500;
                // RTRNewtonSolver->Max_Inner_Iter = 500;
                // ROPTLIB::PARAMSMAP solverParams = {std::pair<std::string, double>("Max_Inner_Iter", 10)};
                // RTRNewtonSolver->SetParams(solverParams);
                RTRNewtonSolver->CheckParams();

                ROFL_VAR1(iterAdmm)
                // startX.Print("startx in ADMM loop");

                // % Solve.
                // [x, xcost, info, options] = trustregions(problem);
                RTRNewtonSolver->Run();
                // Numerically check gradient consistency (optional).
                auto Xopt = RTRNewtonSolver->GetXopt();
                auto XoptCost = RTRNewtonSolver->Getfinalfun();
                costOut = XoptCost;

                // ROFL_VAR1("")
                // Prob.CheckGradHessian(Xopt);

                // std::cout << "Prob.GetUseGrad() " << Prob.GetUseGrad() << std::endl;
                // std::cout << "Prob.GetUseHess() " << Prob.GetUseHess() << std::endl;
                // std::cout << "Prob.GetNumGradHess() " << Prob.GetNumGradHess() << std::endl;

                // params.z = max(ones(size(lambdas_manopt_out)), lambdas_manopt_out);
                // params.y = params.y + params.mu *(params.z-lambdas_manopt_out);

                SomUtils::MatD LambdasManoptOutEig = SomUtils::MatD::Zero(e, 1);
                SomUtils::MatD XoptEig = SomUtils::MatD::Zero(d * d * n + d * n + e, 1);
                ProbAdmm.RoptToEig(Xopt, XoptEig);
                ProbAdmm.getScales(XoptEig, LambdasManoptOutEig);
                // SomUtils::MatD tmp = LambdasManoptOutEig - (yAdmm / muAdmm); // as in Overleaf, but this is not what Matlab code does
                SomUtils::MatD tmp = LambdasManoptOutEig; // as in Matlab
                zAdmm = tmp.cwiseMax(1.0);
                yAdmm += muAdmm * (zAdmm - LambdasManoptOutEig);

                ROFL_VAR1("ADMM updates")
                ROFL_VAR1(zAdmm.transpose())
                ROFL_VAR1(yAdmm.transpose())

                Xopt.CopyTo(startXlocal);

                // Outputs
                Xopt.Print("Xopt");
                std::cout << "XoptCost " << XoptCost << std::endl; // x cost

                delete RTRNewtonSolver;
                // end of lsomRTR()

                // [params.mu, r_k, s_k] = update_lsom_penalty_param(params.mu, x_k, z_k, z_prev);
                // params_mu_next = params.mu;
                // disp(params.mu)
                ProbAdmm.updateLsomPenaltyParam(ProbAdmm.muAdmm_, LambdasManoptOutEig, zAdmm, ProbAdmm.zAdmm_, muAdmm, rK, sK);
                ROFL_VAR1(muAdmm)

                // disp(norm(s_k))
                // disp(norm(r_k))
                ROFL_VAR2(rK.norm(), sK.norm())

                // if size (X_manopt_out.R, 1) == size(X_manopt_out.R, 2)
                //     disp("multidet(X_manopt_out.R)")
                //     disp(multidet(X_manopt_out.R))

                SomUtils::VecMatD RadmmOut(n, SomUtils::MatD::Zero(d, d));
                std::vector<double> rotDetsOk(n, -1.0);
                ProbAdmm.getRotations(XoptEig, RadmmOut);
                SomUtils::multidet(RadmmOut, rotDetsOk);

                // admm_stopping_condition_reached = check_admm_stopping_condition(x_k, y_k, z_k, r_k, s_k, num_edges, 1e-8, 1e-8);
                admmStoppingConditionReached = ProbAdmm.checkAdmmStoppingCondition(rK, sK, ProbAdmm.tolAdmmPrimal_, ProbAdmm.tolAdmmDual_);
                ROFL_VAR1(admmStoppingConditionReached)

                iterAdmm++;
            }
        }
        else
        {
            rsActuallyUseful = true;

            // // RS
            int r0 = d + 1;

            integer numoftypes = 3; // 2 i.e. (3D) Stiefel + Euclidean
            integer numofmani1 = n; // num of Stiefel manifolds
            integer numofmani2 = 1;
            integer numofmani3 = 1;

            ROPTLIB::Stiefel mani1(d, d);
            mani1.ChooseParamsSet2();
            ROPTLIB::Euclidean mani2(d, n);
            ROPTLIB::Euclidean mani3(e);
            ROPTLIB::ProductManifold ProdManiLsom(numoftypes,
                                                  &mani1, numofmani1, &mani2, numofmani2, &mani3, numofmani3);

            SomUtils::MatD XoptEigVec(SomUtils::MatD::Zero(d * d * n + d * n + e, 1));
            Prob.RoptToEig(startXlocal, XoptEigVec);
            // ROFL_VAR3(XoptEigVec.transpose(), XoptEigVec.rows(), e)

            double costLast = std::numeric_limits<double>::infinity();
            auto ProbPrev = Prob;
            // int staircaseStepIdx;
            SomUtils::VecMatD RmanoptOutEig(n, SomUtils::MatD::Zero(d, d));
            SomUtils::MatD TmanoptOutEig(SomUtils::MatD::Zero(d, n));
            SomUtils::MatD LambdaManoptOutEig(SomUtils::MatD::Zero(e, 1));

            for (staircaseStepIdx = r0; staircaseStepIdx <= d * d * n + 1; ++staircaseStepIdx)
            {
                ROFL_VAR1(staircaseStepIdx)
                ROFL_VAR1(costLast)

                SomUtils::VecMatD R(n, SomUtils::MatD::Zero(staircaseStepIdx - 1, d));
                SomUtils::MatD T(SomUtils::MatD::Zero(staircaseStepIdx - 1, n));
                SomUtils::MatD Lambdas(SomUtils::MatD::Zero(e, 1));
                {
                    SomUtils::SomSize somSzScope(staircaseStepIdx - 1, d, n); // TODO: improve getRotations() and getTranslations() and avoid local scope
                    ROPTLIB::LsomProblem ProbScope(somSzScope, Prob.tijs_, Prob.edges_);

                    // auto XoptLocal = XoptNext;
                    ROFL_VAR1("Calling getRotations()")
                    ProbScope.getRotations(XoptEigVec, R);
                    ProbScope.getTranslations(XoptEigVec, T);
                    ProbScope.getScales(XoptEigVec, Lambdas);

                    int iterAdmm = 0;

                    SomUtils::MatD rK(SomUtils::MatD::Zero(e, 1));
                    SomUtils::MatD sK(SomUtils::MatD::Zero(e, 1));

                    // ProbScope.EigToRopt(XoptEigVec, RmanoptOutEig, TmanoptOutEig, LambdaManoptOutEig);

                    // while iter_admm < 100 && ~admm_stopping_condition_reached
                    while (iterAdmm < Prob.maxIterAdmm_ && !admmStoppingConditionReached)
                    {
                        ROFL_VAR2(iterAdmm, Prob.maxIterAdmm_)

                        auto ProbAdmm = Prob;
                        ProbAdmm.costCurr_ = costOut; // TODO: use setter instead of direct access
                        ProbAdmm.setZAdmm(zAdmm);
                        ProbAdmm.setYAdmm(yAdmm);
                        ProbAdmm.setMuAdmm(muAdmm);

                        // lsomRTR(nrs, d, N, problem_data, params, transf_initguess_struct, lambdas_initguess)
                        ROPTLIB::RTRNewton *RTRNewtonSolver = new ROPTLIB::RTRNewton(&ProbAdmm, &startXlocal); // USE INITGUESS HERE!
                        RTRNewtonSolver->Verbose = ROPTLIB::ITERRESULT;
                        // RTRNewtonSolver->Max_Iteration = 500;
                        // RTRNewtonSolver->Max_Inner_Iter = 500;
                        // ROPTLIB::PARAMSMAP solverParams = {std::pair<std::string, double>("Max_Inner_Iter", 10)};
                        // RTRNewtonSolver->SetParams(solverParams);
                        RTRNewtonSolver->CheckParams();

                        ROFL_VAR1(iterAdmm)
                        // startX.Print("startx in ADMM loop");

                        // % Solve.
                        // [x, xcost, info, options] = trustregions(problem);
                        RTRNewtonSolver->Run();
                        // Numerically check gradient consistency (optional).
                        auto Xopt = RTRNewtonSolver->GetXopt();
                        auto XoptCost = RTRNewtonSolver->Getfinalfun();
                        costOut = XoptCost;

                        // ROFL_VAR1("")
                        // Prob.CheckGradHessian(Xopt);

                        // std::cout << "Prob.GetUseGrad() " << Prob.GetUseGrad() << std::endl;
                        // std::cout << "Prob.GetUseHess() " << Prob.GetUseHess() << std::endl;
                        // std::cout << "Prob.GetNumGradHess() " << Prob.GetNumGradHess() << std::endl;

                        // params.z = max(ones(size(lambdas_manopt_out)), lambdas_manopt_out);
                        // params.y = params.y + params.mu *(params.z-lambdas_manopt_out);

                        SomUtils::MatD LambdasManoptOutEig = SomUtils::MatD::Zero(e, 1);
                        SomUtils::MatD XoptEig = SomUtils::MatD::Zero(d * d * n + d * n + e, 1);
                        ProbAdmm.RoptToEig(Xopt, XoptEig);
                        ProbAdmm.getScales(XoptEig, LambdasManoptOutEig);
                        // SomUtils::MatD tmp = LambdasManoptOutEig - (yAdmm / muAdmm); // as in Overleaf, but this is not what Matlab code does
                        SomUtils::MatD tmp = LambdasManoptOutEig; // as in Matlab
                        zAdmm = tmp.cwiseMax(1.0);
                        yAdmm += muAdmm * (zAdmm - LambdasManoptOutEig);

                        ROFL_VAR1("ADMM updates")
                        ROFL_VAR1(zAdmm.transpose())
                        ROFL_VAR1(yAdmm.transpose())

                        Xopt.CopyTo(startXlocal);

                        // Outputs
                        Xopt.Print("Xopt");
                        std::cout << "XoptCost " << XoptCost << std::endl; // x cost

                        delete RTRNewtonSolver;
                        // end of lsomRTR()

                        // [params.mu, r_k, s_k] = update_lsom_penalty_param(params.mu, x_k, z_k, z_prev);
                        // params_mu_next = params.mu;
                        // disp(params.mu)
                        ProbAdmm.updateLsomPenaltyParam(ProbAdmm.muAdmm_, LambdasManoptOutEig, zAdmm, ProbAdmm.zAdmm_, muAdmm, rK, sK);
                        ROFL_VAR1(muAdmm)

                        // disp(norm(s_k))
                        // disp(norm(r_k))
                        ROFL_VAR2(rK.norm(), sK.norm())

                        // if size (X_manopt_out.R, 1) == size(X_manopt_out.R, 2)
                        //     disp("multidet(X_manopt_out.R)")
                        //     disp(multidet(X_manopt_out.R))

                        SomUtils::VecMatD RadmmOut(n, SomUtils::MatD::Zero(d, d));
                        std::vector<double> rotDetsOk(n, -1.0);
                        ProbAdmm.getRotations(XoptEig, RadmmOut);
                        SomUtils::multidet(RadmmOut, rotDetsOk);

                        // admm_stopping_condition_reached = check_admm_stopping_condition(x_k, y_k, z_k, r_k, s_k, num_edges, 1e-8, 1e-8);
                        admmStoppingConditionReached = ProbAdmm.checkAdmmStoppingCondition(rK, sK, ProbAdmm.tolAdmmPrimal_, ProbAdmm.tolAdmmDual_);
                        ROFL_VAR1(admmStoppingConditionReached)

                        iterAdmm++;
                    }
                }

                // SomUtils::VecMatD Rnext(n, SomUtils::MatD::Zero(staircaseStepIdx, d));
                // SomUtils::MatD Tnext(SomUtils::MatD::Zero(staircaseStepIdx, n));
                // SomUtils::MatD LambdasNext(SomUtils::MatD::Zero(e, 1));

                // SomUtils::catZeroRow3dArray(R, Rnext);
                // SomUtils::catZeroRow(T, Tnext);

                SomUtils::SomSize somSzNext(staircaseStepIdx, d, n);
                ROPTLIB::LsomProblem ProbNext(somSzNext, Prob.tijs_, Prob.edges_);

                ProbNext.setGt(rGt, tGt, lambdasGt);
                ProbNext.setRho(Prob.rho_);

                if (Prob.reluScaleCompensation_)
                {
                    ROFL_VAR1(ProbNext.costEigenRelu(ProbNext.Rgt_, ProbNext.Tgt_, ProbNext.LambdasGt_));
                    ROFL_VAR1(ProbNext.costEigenRelu(R, T, Lambdas));
                }
                else
                {
                    ROFL_VAR1(ProbNext.costEigen(ProbNext.Rgt_, ProbNext.Tgt_, ProbNext.LambdasGt_));
                    ROFL_VAR1(ProbNext.costEigen(R, T, Lambdas));
                }

                ROPTLIB::Stiefel mani1next(somSzNext.p_, somSzNext.d_);
                mani1next.ChooseParamsSet2();
                ROPTLIB::Euclidean mani2next(somSzNext.p_, somSzNext.n_);
                ROPTLIB::ProductManifold ProdManiNext(numoftypes,
                                                      &mani1next, numofmani1, &mani2next, numofmani2, &mani3, numofmani3);
                ROPTLIB::Vector Y0;
                SomUtils::VecMatD vR(n, SomUtils::MatD::Zero(somSzNext.p_, somSzNext.d_));
                SomUtils::MatD vT(SomUtils::MatD::Zero(somSzNext.p_, somSzNext.n_));
                SomUtils::MatD vLambdas(SomUtils::MatD::Zero(Prob.numEdges_, 1));
                ProbPrev.setCostCurr(costLast);
                // for (auto &Rm : R)
                //     ROFL_VAR1(ProbPrev.checkIsOnStiefel(Rm))

                double lambda;
                if (Prob.usePIM_)
                {
                    ROFL_VAR1("Calling Prob.lsomEscapeHessianGenprocEigenPIM()")
                    ProbNext.lsomPimHessianGenprocEigen(1e-5, R, T, Lambdas, Y0, lambda, vR, vT, vLambdas); //!! catZeroRows() increase is being done inside
                }
                else
                {
                    ROFL_VAR1("Calling ProbPrev.lsomEscapeHessianGenprocEigen()")
                    ProbPrev.lsomEscapeHessianGenprocEigen(R, T, Lambdas, Y0, lambda, vR, vT, vLambdas);
                }

                if (lambda > -1e-8)
                {
                    ROFL_VAR2(lambda, "R, T eigenvals > 0: exiting stsaircase")
                    // staircaseStepSkipped = 0;
                    RmanoptOutEig = R;
                    TmanoptOutEig = T;
                    LambdaManoptOutEig = Lambdas;

                    costOut = costLast;

                    break;
                }

                Y0.Print("Y0 before costNewStart");
                double costNewStart = ProbNext.f(Y0);
                ROFL_VAR1(costNewStart)

                // Run next step of staircase with found initial guess

                Y0.Print("Y0");

                // Set Prob params
                ProbNext.SetDomain(&ProdManiNext);
                ProbNext.SetUseGrad(true);
                ProbNext.SetUseHess(true);

                ROPTLIB::RTRNewton *RTRNewtonSolverNext = new ROPTLIB::RTRNewton(&ProbNext, &Y0); // USE INITGUESS HERE!
                RTRNewtonSolverNext->Verbose = ROPTLIB::ITERRESULT;
                // RTRNewtonSolverNext->Max_Iteration = 500;
                // RTRNewtonSolverNext->Max_Inner_Iter = 500;
                // ROPTLIB::PARAMSMAP solverParams = {std::pair<std::string, double>("Max_Inner_Iter", 10)};
                // RTRNewtonSolverNext->SetParams(solverParams);
                RTRNewtonSolverNext->CheckParams();

                // Solve.
                // [x, xcost, info, options] = trustregions(problem);
                RTRNewtonSolverNext->Run();
                auto XoptNext = RTRNewtonSolverNext->GetXopt();
                realdp XoptNextCost = RTRNewtonSolverNext->Getfinalfun();
                // Numerically check gradient consistency (optional).
                ProbNext.CheckGradHessian(XoptNext);

                costLast = XoptNextCost;

                ROFL_VAR1(costLast)

                XoptEigVec.resize(staircaseStepIdx * d * n + staircaseStepIdx * n + e, 1);
                ProbNext.RoptToEig(XoptNext, XoptEigVec);
                XoptNext.Print("XoptNext");
                SomUtils::VecMatD XoptNextR(n, SomUtils::MatD::Zero(staircaseStepIdx, d));
                SomUtils::MatD XoptNextT(SomUtils::MatD::Zero(staircaseStepIdx, n));
                SomUtils::MatD XoptNextLambdas(SomUtils::MatD::Zero(e, 1));
                ROFL_VAR1("Calling getRotations()")
                ProbNext.getRotations(XoptEigVec, XoptNextR);
                ProbNext.getTranslations(XoptEigVec, XoptNextT);
                ProbNext.getScales(XoptEigVec, XoptNextLambdas);

                // std::cout << "Prob.GetUseGrad() " << Prob.GetUseGrad() << std::endl;
                // std::cout << "Prob.GetUseHess() " << Prob.GetUseHess() << std::endl;
                // std::cout << "Prob.GetNumGradHess() " << Prob.GetNumGradHess() << std::endl;

                // Outputs
                XoptNext.Print("XoptNext");
                XoptNext.CopyTo(startXlocal);
                std::cout << "XoptNextCost " << XoptNextCost << std::endl; // x cost

                ProbNext.RoptToEig(XoptNext, XoptEigVec);
                ROFL_VAR1(XoptEigVec.transpose())

                ProbPrev = ProbNext;

                // save output
                for (int i = 0; i < n; ++i)
                {
                    RmanoptOutEig[i].resize(staircaseStepIdx, d);
                }
                TmanoptOutEig.resize(staircaseStepIdx, n);

                RmanoptOutEig = XoptNextR;
                TmanoptOutEig = XoptNextT;
                LambdaManoptOutEig = XoptNextLambdas;

                // // Rank stopping condition
                // ROFL_VAR1(staircaseStepIdx)
                // SomUtils::MatD XoutRhSt(SomUtils::MatD::Zero(staircaseStepIdx, d * n));
                // SomUtils::hstack(RmanoptOutEig, XoutRhSt);
                // Eigen::FullPivLU<SomUtils::MatD> luDecomp(XoutRhSt);
                // auto rank = luDecomp.rank();
                // if (rank < staircaseStepIdx)
                // {
                //     staircaseStepIdx++;
                //     ROFL_VAR1("Rank stopping condition reached -> Exiting RS");
                //     break;
                // }

                delete RTRNewtonSolverNext;

                // break; // uncomment this to run only one step of staircase
            }

            // // Recovery procedure

            ROFL_VAR1("Running recovery procedure")

            // back to SE(d)^N

            for (int i = 0; i < n; ++i)
            {
                ROFL_VAR2(i, RmanoptOutEig[i])
            }
            ROFL_VAR1(TmanoptOutEig)
            ROFL_VAR1(LambdaManoptOutEig)

            ROFL_VAR1(staircaseStepIdx) // unused in ProbPrev.recoverySEdN() anyway

            SomUtils::VecMatD Rrecovered(n, SomUtils::MatD::Zero(d, d));
            SomUtils::MatD Trecovered(SomUtils::MatD::Zero(d, n));
            SomUtils::MatD LambdasRecovered(SomUtils::MatD::Zero(e, 1));
            bool recSEDNsuccess = ProbPrev.recoverySEdN(staircaseStepIdx,
                                                        RmanoptOutEig, TmanoptOutEig, LambdaManoptOutEig,
                                                        Rrecovered, Trecovered, LambdasRecovered);

            if (!recSEDNsuccess)
            {
                ROFL_VAR1("Recovery procedure failed")
                // SomUtils::MatD XbackupEig(SomUtils::MatD::Zero(d * d * n + d * n + e, 1));
                // Prob.RoptToEig(Xbackup, XbackupEig);
                // Prob.getRotations(XbackupEig, Rrecovered);
                // Prob.getTranslations(XbackupEig, Trecovered);
                // Prob.getScales(XbackupEig, LambdasRecovered);
                Xbackup.CopyTo(startXlocal);
                rsActuallyUseful = false; // TODO: set this flag to true if recovery procedure is actually useful
            }

            ROFL_VAR1("Printing R, T, lambdas after recovery")
            for (auto &m : Rrecovered)
                ROFL_VAR1(m)
            ROFL_VAR1(Trecovered)
            ROFL_VAR1(LambdasRecovered)

            ROFL_VAR1(recSEDNsuccess)
            rsSuccess = recSEDNsuccess;
        }

        ROFL_VAR1("End of ADMM")

        //
        std::vector<double> Rdets(n);
        SomUtils::VecMatD Rrecovered(n, SomUtils::MatD::Zero(d, d));
        SomUtils::MatD Trecovered(SomUtils::MatD::Zero(d, n));
        SomUtils::MatD LambdasRecovered(SomUtils::MatD::Zero(e, 1));
        SomUtils::VecMatD RmanoptOutEig(n, SomUtils::MatD::Zero(staircaseStepIdx, d));
        SomUtils::MatD TmanoptOutEig(SomUtils::MatD::Zero(staircaseStepIdx, n));
        SomUtils::MatD LambdaManoptOutEig(SomUtils::MatD::Zero(e, 1));
        SomUtils::MatD XoptEig(SomUtils::MatD::Zero(d * d * n + d * n + e, 1));
        Prob.RoptToEig(startXlocal, XoptEig);
        Prob.getRotations(XoptEig, Rrecovered);
        Prob.getTranslations(XoptEig, Trecovered);
        Prob.getScales(XoptEig, LambdasRecovered);
        SomUtils::multidet(Rrecovered, Rdets);

        rotDetsOk = true;
        for (int i = 0; i < n; ++i)
        {
            // ROFL_VAR2(i, Rdets[i])
            if (!SomUtils::isEqualDoubles(fabs(Rdets[i]), 1))
            {
                rotDetsOk = false;
                ROFL_VAR1("Rotation determinant condition failed")
                break;
            }
        }

        lambdasAcceptable = true;
        for (int i = 0; i < e; ++i)
        {
            if (LambdasRecovered(i, 0) < 1)
            {
                lambdasAcceptable = false;
                ROFL_VAR1("Lambda negativity condition failed")
                break;
            }
        }

        // globalize

        Rout.resize(n, SomUtils::MatD::Zero(d, d));
        Tout.resize(d, n);
        Tout.setZero();
        lambdasOut.resize(e, 1);
        lambdasOut.setZero();

        bool globalRecoverySuccess = true; // TODO: implement globalization procedure and set this flag accordingly

        if (Prob.performGlobalization_)
        {
            ROFL_VAR1("Running globalization procedure")
            bool globalRecoverySuccess = Prob.globalize(src, Rrecovered, Trecovered, LambdasRecovered,
                                                        Rout, Tout, lambdasOut);
            if (!globalRecoverySuccess)
            {
                ROFL_VAR1("Globalization procedure failed")
                Rout = Rrecovered;
                Tout = Trecovered;
                lambdasOut = LambdasRecovered;

                rsActuallyUseful = false; // TODO: set this flag to true if globalization procedure is actually useful

                SomUtils::MatD XbackupEig(SomUtils::MatD::Zero(d * d * n + d * n + e, 1));
                SomUtils::VecMatD Rbackup(n, SomUtils::MatD::Zero(d, d));
                SomUtils::MatD Tbackup(SomUtils::MatD::Zero(d, n));
                SomUtils::MatD LambdasBackup(SomUtils::MatD::Zero(e, 1));
                Prob.RoptToEig(Xbackup, XbackupEig);
                Prob.getRotations(XbackupEig, Rbackup);
                Prob.getTranslations(XbackupEig, Tbackup);
                Prob.getScales(XbackupEig, LambdasBackup);
                bool globalRecoverySuccess = Prob.globalize(src, Rbackup, Tbackup, LambdasBackup,
                                                            Rout, Tout, lambdasOut);
            }
        }
        else
        {
            ROFL_VAR1("Skipping globalization procedure")
            Rout = Rrecovered;
            Tout = Trecovered;
            lambdasOut = LambdasRecovered;
        }

        for (int i = 0; i < n; ++i)
        {
            ROFL_VAR2(i, Rout[i])
            if (!SomUtils::isEqualDoubles(((Rout[i] - rGt[i]).cwiseAbs().maxCoeff()), 0.0))
            {
                rotDetsOk = false;
                ROFL_VAR1("Rotations changed")
                break;
            }
        }
        ROFL_VAR1(Tout)

        if (!SomUtils::isEqualDoubles(((Tout - tGt).cwiseAbs().maxCoeff()), 0.0))
        {
            rotDetsOk = false;
            ROFL_VAR1("Translations changed")
        }

        ROFL_VAR1(lambdasOut)
        if (!SomUtils::isEqualDoubles(((lambdasOut - lambdasGt).cwiseAbs().maxCoeff()), 0.0))
        {
            rotDetsOk = false;
            ROFL_VAR1("Lambdas changed")
        }

        ROFL_VAR1(globalRecoverySuccess)

        return costOut;
    }

} // end of namespace ROPTLIB
