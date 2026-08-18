#include "problems_Ssom.h"

namespace ROPTLIB
{

    double runSsomICP(ROPTLIB::SsomProblem &Prob,
                      const ROPTLIB::Vector &startX,
                      int src,
                      SomUtils::VecMatD &Rout,
                      SomUtils::MatD &Tout,
                      SomUtils::MatD &LambdasOut,
                      double stopThr)
    {
        int d = Prob.sz_.d_;
        int n = Prob.sz_.n_;
        int e = Prob.edges_.rows();

        ROFL_VAR1("Start of runRsomRS()")
        ROFL_VAR1(Prob.costEigen(Prob.Rgt_, Prob.Tgt_, Prob.LambdasGt_));

        // output the parameters of the manifold of domain
        ROPTLIB::RTRNewton *RTRNewtonSolver = new ROPTLIB::RTRNewton(&Prob, &startX); // USE INITGUESS HERE!
        RTRNewtonSolver->Verbose = ROPTLIB::ITERRESULT;
        // RTRNewtonSolver->Max_Iteration = 500;
        // RTRNewtonSolver->Max_Inner_Iter = 500;
        // ROPTLIB::PARAMSMAP solverParams = {std::pair<std::string, double>("Max_Inner_Iter", 10)};
        // RTRNewtonSolver->SetParams(solverParams);
        RTRNewtonSolver->CheckParams();

        // % Solve.
        // [x, xcost, info, options] = trustregions(problem);
        RTRNewtonSolver->Run();
        // Numerically check gradient consistency (optional).
        auto Xopt = RTRNewtonSolver->GetXopt();
        auto XoptCost = RTRNewtonSolver->Getfinalfun();

        // Prob.CheckGradHessian(Xopt);

        // std::cout << "Prob.GetUseGrad() " << Prob.GetUseGrad() << std::endl;
        // std::cout << "Prob.GetUseHess() " << Prob.GetUseHess() << std::endl;
        // std::cout << "Prob.GetNumGradHess() " << Prob.GetNumGradHess() << std::endl;

        // Outputs
        Xopt.Print("Xopt");
        std::cout << "XoptCost " << XoptCost << std::endl; // xopt cost

        delete RTRNewtonSolver;

        // ICP

        SomUtils::MatD XoptEigVec(SomUtils::MatD::Zero(d * d * n + d * n + e, 1));
        Prob.RoptToEig(Xopt, XoptEigVec);
        ROFL_VAR1(XoptEigVec.transpose())

        double costLast = XoptCost;
        auto rGt = Prob.Rgt_;
        auto tGt = Prob.Tgt_;
        auto lambdasGt = Prob.LambdasGt_;


        SomUtils::MatD startxEigVec(SomUtils::MatD::Zero(d * d * n + d * n + e, 1));
        Prob.RoptToEig(startX, startxEigVec);
        double diffPrevCurr = (startxEigVec - XoptEigVec).norm();

        auto XoptPrev = Xopt;
        SomUtils::MatD XoptPrevEigVec(SomUtils::MatD::Zero(d * d * n + d * n + e, 1));

        for (int i = 0; i < Prob.icpMaxIterations_; ++i)
        {
            ROFL_VAR1(i)
            if (diffPrevCurr < stopThr) // TODO: make 1e-3 a parameter
            {
                ROFL_VAR1("diffPrevCurr < stopThr -> stopping condition reached")
                Prob.RoptToEig(XoptPrev, XoptPrevEigVec);
                break;
            }

            ROPTLIB::RTRNewton *RTRNewtonSolverIter = new ROPTLIB::RTRNewton(&Prob, &XoptPrev); // USE INITGUESS HERE!
            RTRNewtonSolverIter->Verbose = ROPTLIB::ITERRESULT;
            RTRNewtonSolverIter->Run();
            auto XoptIter = RTRNewtonSolverIter->GetXopt();
            auto XoptIterCost = RTRNewtonSolverIter->Getfinalfun();

            // TODO: use
            // Prob.ssomStep3(rGt, tGt, LambdasOut); 
            // instead of the full RTR solution

            // Outputs
            XoptIter.Print("XoptIter");
            std::cout << "XoptIterCost " << XoptIterCost << std::endl; // x cost

            delete RTRNewtonSolverIter;

            SomUtils::MatD XoptIterEigVec(SomUtils::MatD::Zero(d * d * n + d * n + e, 1));
            Prob.RoptToEig(XoptIter, XoptIterEigVec);
            Prob.RoptToEig(XoptPrev, XoptPrevEigVec);
            diffPrevCurr = (XoptPrevEigVec - XoptIterEigVec).norm();

            XoptPrev = XoptIter;
        }

        // globalize

        SomUtils::VecMatD Rlocal(n, SomUtils::MatD::Zero(d, d));
        SomUtils::MatD Tlocal(SomUtils::MatD::Zero(d, n));
        SomUtils::MatD LambdasLocal(SomUtils::MatD::Zero(e, 1));
        Prob.RoptToEig(XoptPrev, XoptPrevEigVec); // in case all iterations are performed

        Prob.getRotations(XoptPrevEigVec, Rlocal); // TODO: improve getRotations() and getTranslations() and stop using them as class methods
        Prob.getTranslations(XoptPrevEigVec, Tlocal);
        Prob.getScales(XoptPrevEigVec, LambdasLocal);

        Rout.resize(n, SomUtils::MatD::Zero(d, d));
        Tout.resize(d, n);

        ROFL_VAR1(Prob.costEigen(Rlocal, Tlocal, LambdasLocal));

        ROFL_VAR1("Solving relative gauge ambiguity")

        SomUtils::VecMatD Rrecovered(n, SomUtils::MatD::Zero(d, d));
        SomUtils::MatD Trecovered(SomUtils::MatD::Zero(d, n));
        SomUtils::MatD LambdasRecovered(SomUtils::MatD::Zero(e, 1));
        bool recSEdn = Prob.recoverySEdN(d + 1,
                                         Rlocal, Tlocal, LambdasLocal,
                                         Rrecovered, Trecovered, LambdasRecovered);
        ROFL_VAR1(recSEdn)
        ROFL_VAR1(Prob.costEigen(Rrecovered, Trecovered, LambdasRecovered));

        ROFL_VAR1("Running globalization procedure")

        bool globalRecoverySuccess = Prob.globalize(src, Rrecovered, Trecovered, LambdasRecovered,
                                                    Rout, Tout, LambdasOut);
        ROFL_VAR1(globalRecoverySuccess)
        ROFL_VAR1(Prob.costEigen(Rout, Tout, LambdasOut));

        return costLast;
    }

    void SsomProblem::ssomStep3(const SomUtils::VecMatD &Rin,
                                const SomUtils::MatD &Tin,
                                SomUtils::MatD &LambdasOut,
                                double stopThr)
    {
        // num_edges = size(edges, 1);
        // lambdas_out = zeros(num_edges, 1);
        int numEdges = edges_.rows();
        // for ee = 1 : num_edges

        //      tij = tijs( :, ee);

        //      ii = edges(ee, 1);
        //      jj = edges(ee, 2);
        //      Ri = R( :, :, ii);
        //      % Rj = R( :, :, jj);
        //      Ti = T( :, ii);
        //      Tj = T( :, jj);

        //      a = tij' * tij; b1 = tij ' * Ri' * Ti;
        //      b2 = -tij ' * Ri' * Tj;
        //      b3 = Ti' * Ri * tij; b4 = -Tj' * Ri * tij; b = b1 + b2 + b3 + b4;
        //      c1 = Ti' * Ti; c2 = -Ti' * Tj; c3 = -Tj' * Ti; % OBS: should be equal to c2 c4 = Tj' * Tj; c = c1 + c2 + c3 + c4;

        //      x = solve_quadratic_lambdas(a, b, c);

        //      lambdas_out(ee) = x;

        for (int e = 0; e < numEdges; ++e)
        {
            auto tij = tijs_.col(e);
            int ii = edges_(e, 0) - 1;
            int jj = edges_(e, 1) - 1;
            auto Ri = Rin[ii];
            auto Ti = Tin.col(ii);
            auto Tj = Tin.col(jj);

            double a = tij.transpose() * tij;
            double b1 = tij.transpose() * Ri.transpose() * Ti;
            double b2 = -tij.transpose() * Ri.transpose() * Tj;
            double b3 = Ti.transpose() * Ri * tij;
            double b4 = -Tj.transpose() * Ri * tij;
            double b = b1 + b2 + b3 + b4;

            double c1 = Ti.transpose() * Ti;
            double c2 = -Ti.transpose() * Tj;
            double c3 = -Tj.transpose() * Ti;
            double c4 = Tj.transpose() * Tj;
            double c = c1 + c2 + c3 + c4;

            LambdasOut(e) = solveQuadraticLambdas(a, b, c);
        }
    }

    double SsomProblem::solveQuadraticLambdas(const double a, const double b, const double c)
    {
        // D = b.^2 - 4.*a.*c;                    % discriminant

        double D = b * b - 4.0 * a * c;

        // if D < 0
        //     x = 1;
        //     return;
        if (D < 0)
        {
            return 1.0;
        }
        else
        {
            double sqrtD = sqrt(D);
            double x1 = (-b + sqrtD) / (2.0 * a);
            double x2 = (-b - sqrtD) / (2.0 * a);
            if (x1 < 1.0 && x2 >= 1.0)
            {
                return x2;
            }
            else if (x2 < 1.0 && x1 >= 1.0)
            {
                return x1;
            }
            else if (x1 >= 1.0 && x2 >= 1.0)
            {
                return std::max(x1, x2);
            }
            else
            {
                // which case are we in?
                return 1.0;
            }
        }
        // else
        //     sqrtD = sqrt(D);
        //     x1 = (-b + sqrtD) / (2.*a);
        //     x2 = (-b - sqrtD) / (2.*a);
        //     if x1 < 1 && x2 >=1
        //         x = x2;
        //     elseif  x2 < 1 && x1 >=1
        //         x = x1;
        //     elseif x1 >= 1 && x2 >= 1
        //         x = max(x1, x2);
        //     else
        //         %which case are we in?
        //         x = 1;
    }

} // end of namespace ROPTLIB