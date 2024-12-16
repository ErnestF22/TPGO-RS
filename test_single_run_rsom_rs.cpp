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

    for (int staircaseStepIdx = r0; staircaseStepIdx <= d * n + 1; ++staircaseStepIdx)
    {
        ROFL_VAR1(staircaseStepIdx)

        SomUtils::VecMatD R(n, SomUtils::MatD::Zero(staircaseStepIdx - 1, d));
        SomUtils::MatD T(SomUtils::MatD::Zero(staircaseStepIdx - 1, n));
        Prob.getRotations(XoptEigVec, R);
        Prob.getTranslations(XoptEigVec, T);

        SomUtils::VecMatD Rnext(n, SomUtils::MatD::Zero(staircaseStepIdx, d));
        SomUtils::MatD Tnext(SomUtils::MatD::Zero(staircaseStepIdx, n));

        Prob.catZeroRow3dArray(R, Rnext);
        Prob.catZeroRow(T, Tnext);

        auto somSzNext = Prob.sz_;
        somSzNext.p_++;

        ROPTLIB::SampleSomProblem ProbNext(somSzNext, Prob.Tijs_, Prob.edges_);

        ROPTLIB::Stiefel mani1next(somSzNext.p_, somSzNext.d_);
        mani1next.ChooseParamsSet2();
        ROPTLIB::Euclidean mani2next(somSzNext.p_, somSzNext.n_);
        ROPTLIB::ProductManifold ProdManiNext(numoftypes, &mani1next, numofmani1, &mani2next, numofmani2);
        ROPTLIB::Vector Y0 = ProdManiNext.RandominManifold();
        ProbNext.rsomPimHessianGenproc(1e-4, Rnext, Tnext, Y0);

        // Run next step of staircase with found initial guess

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
        // Numerically check gradient consistency (optional).
        auto XoptNext = RTRNewtonSolverNext->GetXopt();
        auto XoptNextCost = RTRNewtonSolverNext->Getfinalfun();

        ProbNext.CheckGradHessian(XoptNext);

        // std::cout << "Prob.GetUseGrad() " << Prob.GetUseGrad() << std::endl;
        // std::cout << "Prob.GetUseHess() " << Prob.GetUseHess() << std::endl;
        // std::cout << "Prob.GetNumGradHess() " << Prob.GetNumGradHess() << std::endl;

        // Outputs
        XoptNext.Print("XoptNext");
        std::cout << "XoptNextCost " << XoptNextCost << std::endl; // x cost

        XoptEigVec.resize(staircaseStepIdx * d * n + staircaseStepIdx * n, 1);
        Prob.RoptToEig(XoptNext, XoptEigVec);
        ROFL_VAR1(XoptEigVec.transpose())

        delete RTRNewtonSolverNext;

        break;
    }

    // Implement recovery procedure
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
    // Read from csv
    ROPTLIB::Vector xGt = ProdMani.RandominManifold();
    SomUtils::readCsvInitguess("../matlab/data/cpp_testdata/tdata_n10_mindeg3/X_gt.csv", xGt);
    xGt.Print("xGt");
    ROPTLIB::SampleSomProblem Prob(somSzD, Tijs, edges);
    // Set the domain of the problem to be the product of Stiefel manifolds
    Prob.SetDomain(&ProdMani);

    // Set Prob params
    Prob.SetUseGrad(true);
    Prob.SetUseHess(true);

    // ROPT to Eig (GT)
    SomUtils::MatD XgtVecEig(SomUtils::MatD::Zero(d * d * n + d * n, 1));
    SomUtils::VecMatD RgtEig(n, SomUtils::MatD::Zero(d, d));
    SomUtils::MatD TgtEig(SomUtils::MatD::Zero(d, n));
    Prob.RoptToEig(xGt, XgtVecEig);
    Prob.getRotations(XgtVecEig, RgtEig);
    Prob.getTranslations(XgtVecEig, TgtEig);
    Prob.setGt(RgtEig, TgtEig);

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