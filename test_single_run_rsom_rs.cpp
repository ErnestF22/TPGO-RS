#include <iostream>
#include <fstream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <filesystem>
#include <boost/lexical_cast.hpp>

#include "thirdparty/roptlib/problems_SampleSom.h"

void runRsomRS(ROPTLIB::SampleSomProblem &Prob, const ROPTLIB::Vector &startX)
{
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

    Prob.CheckGradHessian(Xopt);

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
    auto ProbPrev = Prob;
    int staircaseStepIdx;
    SomUtils::VecMatD RmanoptOutEig(n, SomUtils::MatD::Zero(d, d));
    SomUtils::MatD TmanoptOutEig(SomUtils::MatD::Zero(d, n));

    auto rGt = Prob.Rgt_;
    auto tGt = Prob.Tgt_;

    int staircaseStepSkipped;

    for (staircaseStepIdx = r0; staircaseStepIdx <= d * n + 1; ++staircaseStepIdx)
    {
        staircaseStepSkipped = 1;
        ROFL_VAR1(staircaseStepIdx)
        ROFL_VAR1(costLast)

        SomUtils::VecMatD R(n, SomUtils::MatD::Zero(staircaseStepIdx - 1, d));
        SomUtils::MatD T(SomUtils::MatD::Zero(staircaseStepIdx - 1, n));
        Prob.getRotations(XoptEigVec, R);
        Prob.getTranslations(XoptEigVec, T);

        // SomUtils::VecMatD Rnext(n, SomUtils::MatD::Zero(staircaseStepIdx, d));
        // SomUtils::MatD Tnext(SomUtils::MatD::Zero(staircaseStepIdx, n));

        // Prob.catZeroRow3dArray(R, Rnext);
        // Prob.catZeroRow(T, Tnext);

        SomUtils::SomSize somSzNext(staircaseStepIdx, d, n);
        ROPTLIB::SampleSomProblem ProbNext(somSzNext, Prob.Tijs_, Prob.edges_);

        ProbNext.setGt(rGt, tGt);

        ROPTLIB::Stiefel mani1next(somSzNext.p_, somSzNext.d_);
        mani1next.ChooseParamsSet2();
        ROPTLIB::Euclidean mani2next(somSzNext.p_, somSzNext.n_);
        ROPTLIB::ProductManifold ProdManiNext(numoftypes, &mani1next, numofmani1, &mani2next, numofmani2);
        ROPTLIB::Vector Y0 = ProdManiNext.RandominManifold();
        double lambda;
        SomUtils::VecMatD vLambdaR(n, SomUtils::MatD::Zero(somSzNext.p_, somSzNext.d_));
        SomUtils::MatD vLambdaT(SomUtils::MatD::Zero(somSzNext.p_, somSzNext.n_));
        ProbPrev.setCostCurr(costLast);
        ProbPrev.rsomPimHessianGenproc(1e-4, R, T, Y0, lambda, vLambdaR, vLambdaT);

        if (lambda < 0)
        {
            ROFL_VAR1("R, T eigenvals > 0: exiting staircase")
            staircaseStepSkipped = 0;
            break;
        }

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

        // % Solve.
        // [x, xcost, info, options] = trustregions(problem);
        RTRNewtonSolverNext->Run();
        auto XoptNext = RTRNewtonSolverNext->GetXopt();
        auto XoptNextCost = RTRNewtonSolverNext->Getfinalfun();
        // Numerically check gradient consistency (optional).
        ProbNext.CheckGradHessian(XoptNext);

        costLast = XoptNextCost;

        ROFL_VAR1(costLast)

        XoptEigVec.resize(staircaseStepIdx * d * n + staircaseStepIdx * n, 1);
        ProbNext.RoptToEig(XoptNext, XoptEigVec);
        XoptNext.Print("XoptNext");
        SomUtils::VecMatD XoptNextR(n, SomUtils::MatD::Zero(staircaseStepIdx, d));
        SomUtils::MatD XoptNextT(SomUtils::MatD::Zero(staircaseStepIdx, n));
        ProbNext.getRotations(XoptEigVec, XoptNextR);
        ProbNext.getTranslations(XoptEigVec, XoptNextT);

        // std::cout << "Prob.GetUseGrad() " << Prob.GetUseGrad() << std::endl;
        // std::cout << "Prob.GetUseHess() " << Prob.GetUseHess() << std::endl;
        // std::cout << "Prob.GetNumGradHess() " << Prob.GetNumGradHess() << std::endl;

        // Outputs
        XoptNext.Print("XoptNext");
        std::cout << "XoptNextCost " << XoptNextCost << std::endl; // x cost

        ProbNext.RoptToEig(XoptNext, XoptEigVec);
        ROFL_VAR1(XoptEigVec.transpose())

        delete RTRNewtonSolverNext;

        ProbPrev = ProbNext;

        // save output
        for (int i = 0; i < n; ++i)
        {
            RmanoptOutEig[i].resize(staircaseStepIdx, d);
        }
        TmanoptOutEig.resize(staircaseStepIdx, n);

        RmanoptOutEig = XoptNextR;
        TmanoptOutEig = XoptNextT;

        // Rank stopping condition
        ROFL_VAR1(staircaseStepIdx)
        SomUtils::MatD XoutRhSt(SomUtils::MatD::Zero(staircaseStepIdx, d * n));
        ProbNext.hstack(RmanoptOutEig, XoutRhSt);
        Eigen::FullPivLU<SomUtils::MatD> luDecomp(XoutRhSt);
        auto rank = luDecomp.rank();
        if (rank < staircaseStepIdx)
        {
            staircaseStepIdx++;
            ROFL_VAR1("Rank stopping condition reached -> Exiting RS");
            break;
        }

        break; // TODO: For now, just 1 RS step allowed -> remove it later after fixing linesearch
    }

    // Recovery procedure

    ROFL_VAR1("Running recovery procedure")

    // back to SE(d)^N
    SomUtils::VecMatD Rrecovered(n, SomUtils::MatD::Zero(d, d));
    SomUtils::MatD Trecovered(SomUtils::MatD::Zero(d, n));
    bool recSEDNsuccess = ProbPrev.recoverySEdN(staircaseStepIdx - staircaseStepSkipped + 1,
                                                RmanoptOutEig, TmanoptOutEig,
                                                Rrecovered, Trecovered);

    ROFL_VAR1("Printing R, T SE(d)^N")
    for (auto &m : Rrecovered)
        ROFL_VAR1(m)
    ROFL_VAR1(Trecovered)

    ROFL_VAR1(recSEDNsuccess)

    // globalize

    ROFL_VAR1("Running globalization procedure")

    SomUtils::VecMatD Rout(n, SomUtils::MatD::Zero(d, d));
    SomUtils::MatD Tout(SomUtils::MatD::Zero(d, n));
    int src = 0;
    bool globalRecoverySuccess = ProbPrev.globalize(src, Rrecovered, Trecovered,
                                                    Rout, Tout);

    ROFL_VAR1("Printing R, T out")
    for (auto &m : Rout)
        ROFL_VAR1(m)
    ROFL_VAR1(Tout)

    ROFL_VAR1(globalRecoverySuccess)
}

int main(int argc, char **argv)
{
    int d = 3;
    int n;
    SomUtils::readSingleIntCsv("../matlab/data/cpp_testdata/tdata_n10_mindeg3/n.csv", n);

    int numEdges;

    SomUtils::readSingleIntCsv("../matlab/data/cpp_testdata/tdata_n10_mindeg3/e.csv", numEdges);

    ROFL_VAR2(n, numEdges); // n = 10, numEdges = 48

    SomUtils::SomSize somSzD(d, d, n);
    SomUtils::MatD Tijs(d, numEdges);
    Eigen::MatrixXi edges(numEdges, 2);

    SomUtils::readMatlabCsvTijs("../matlab/data/cpp_testdata/tdata_n10_mindeg3/tijs.csv", Tijs, d, numEdges);
    SomUtils::readMatlabCsvEdges("../matlab/data/cpp_testdata/tdata_n10_mindeg3/edges.csv", edges);

    ROFL_VAR1(Tijs)
    ROFL_VAR1(edges)

    // problem d x d x n
    integer numoftypes = 2; // 2 i.e. (3D) Stiefel + Euclidean
    integer numofmani1 = n; // num of Stiefel manifolds
    integer numofmani2 = 1;
    ROPTLIB::Stiefel mani1(d, d);
    mani1.ChooseParamsSet2();
    ROPTLIB::Euclidean mani2(d, n);
    ROPTLIB::ProductManifold ProdMani(numoftypes, &mani1, numofmani1, &mani2, numofmani2);
    ROPTLIB::SampleSomProblem Prob(somSzD, Tijs, edges);

    // Read GT from csv
    ROPTLIB::Vector xGt = ProdMani.RandominManifold();
    SomUtils::readCsvInitguess("../matlab/data/cpp_testdata/tdata_n10_mindeg3/X_gt.csv", xGt);
    xGt.Print("xGt");
    // ROPT to Eig (GT)
    SomUtils::MatD XgtVecEig(SomUtils::MatD::Zero(d * d * n + d * n, 1));
    SomUtils::VecMatD RgtEig(n, SomUtils::MatD::Zero(d, d));
    SomUtils::MatD TgtEig(SomUtils::MatD::Zero(d, n));
    Prob.RoptToEig(xGt, XgtVecEig);
    Prob.getRotations(XgtVecEig, RgtEig);
    Prob.getTranslations(XgtVecEig, TgtEig);
    Prob.setGt(RgtEig, TgtEig);

    // Set the domain of the problem to be the product of Stiefel manifolds
    Prob.SetDomain(&ProdMani);

    // Set Prob params
    Prob.SetUseGrad(true);
    Prob.SetUseHess(true);

    ROFL_VAR1("Printing R, T gt")
    for (auto &m : RgtEig)
        ROFL_VAR1(m)
    ROFL_VAR1(TgtEig)

    // // problem nrs x d x n
    // ROPTLIB::Stiefel mani1nrs(nrs, d);
    // mani1.ChooseParamsSet2();
    // ROPTLIB::Euclidean mani2nrs(nrs, n);
    // ROPTLIB::ProductManifold ProdManiNrs(numoftypes, &mani1nrs, numofmani1, &mani2nrs, numofmani2);
    // // Read from csv
    // ROPTLIB::Vector xManoptOut = ProdManiNrs.RandominManifold();
    // SomUtils::readCsvInitguess("../data/recov_x_manopt_out.csv", xManoptOut);
    // xManoptOut.Print("xManoptOut");

    // Generate startX (random)
    ROPTLIB::Vector startX = ProdMani.RandominManifold();

    // RUN RSOM RS
    runRsomRS(Prob, startX);
}