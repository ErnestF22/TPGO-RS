#include "problems_SampleSom.h"

namespace ROPTLIB
{

    double runRsomICP(SampleSomProblem &Prob, const Vector &startX, int src, SomUtils::VecMatD &Rout, SomUtils::MatD &Tout)
    {
        ROFL_VAR1("Start of runRsomRS()")
        ROFL_VAR1(Prob.costEigen(Prob.Rgt_, Prob.Tgt_));

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
        std::cout << "XoptCost " << XoptCost << std::endl; // x cost

        delete RTRNewtonSolver;

        // RS
        int d = Prob.sz_.d_;
        int n = Prob.sz_.n_;
        int r0 = d + 1;

        integer numoftypes = 2; // 2 i.e. (3D) Stiefel + Euclidean
        integer numofmani1 = n; // num of Stiefel manifolds
        integer numofmani2 = 1;

        SomUtils::MatD XoptEigVec(SomUtils::MatD::Zero(d * d * n + d * n, 1));
        Prob.RoptToEig(Xopt, XoptEigVec);
        ROFL_VAR1(XoptEigVec.transpose())

        double costLast = XoptCost;
        auto rGt = Prob.Rgt_;
        auto tGt = Prob.Tgt_;

        SomUtils::MatD startxEigVec(SomUtils::MatD::Zero(d * d * n + d * n, 1));
        Prob.RoptToEig(startX, startxEigVec);
        double diffPrevCurr = (startxEigVec - XoptEigVec).norm();

        auto XoptPrev = Xopt;
        SomUtils::MatD XoptPrevEigVec(SomUtils::MatD::Zero(d * d * n + d * n, 1));

        for (int i = 0; i < 20; ++i) // TODO: make 20 a parameter
        {
            ROFL_VAR1(i)
            if (diffPrevCurr < 1e-3) // TODO: make 1e-3 a parameter
            {
                Prob.RoptToEig(XoptPrev, XoptPrevEigVec);
                break;
            }

            ROPTLIB::RTRNewton *RTRNewtonSolverIter = new ROPTLIB::RTRNewton(&Prob, &XoptPrev); // USE INITGUESS HERE!
            RTRNewtonSolverIter->Verbose = ROPTLIB::ITERRESULT;
            RTRNewtonSolverIter->Run();
            auto XoptIter = RTRNewtonSolverIter->GetXopt();
            auto XoptIterCost = RTRNewtonSolverIter->Getfinalfun();

            // Outputs
            XoptIter.Print("XoptIter");
            std::cout << "XoptIterCost " << XoptIterCost << std::endl; // x cost

            delete RTRNewtonSolverIter;

            SomUtils::MatD XoptIterEigVec(SomUtils::MatD::Zero(d * d * n + d * n, 1));
            Prob.RoptToEig(XoptIter, XoptIterEigVec);
            Prob.RoptToEig(XoptPrev, XoptPrevEigVec);
            diffPrevCurr = (XoptPrevEigVec - XoptIterEigVec).norm();

            XoptPrev = XoptIter;
        }

        // globalize

        SomUtils::VecMatD Rlocal(n, SomUtils::MatD::Zero(d, d));
        SomUtils::MatD Tlocal(SomUtils::MatD::Zero(d, n));
        Prob.RoptToEig(XoptPrev, XoptPrevEigVec); // in case all iterations are performed

        Prob.getRotations(XoptPrevEigVec, Rlocal); // TODO: improve getRotations() and getTranslations() and stop using them as class methods
        Prob.getTranslations(XoptPrevEigVec, Tlocal);

        ROFL_VAR1("Running globalization procedure")

        Rout.resize(n, SomUtils::MatD::Zero(d, d));
        Tout.resize(d, n);
        bool globalRecoverySuccess = Prob.globalize(src, Rlocal, Tlocal,
                                                    Rout, Tout);
        ROFL_VAR1(globalRecoverySuccess)

        return costLast;
    }

} // end of namespace ROPTLIB