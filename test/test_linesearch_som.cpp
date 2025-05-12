
#include <iostream>
#include <fstream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <filesystem>
#include <boost/lexical_cast.hpp>

#include "thirdparty/roptlib/problems_SampleSom.h"

using namespace ROPTLIB;

void testSomSample(SomUtils::SomSize somSz, SomUtils::MatD &Tijs, Eigen::MatrixXi &edges)
{
    ROFL_VAR3(somSz.p_, somSz.d_, somSz.n_);

    // manifold
    // An example of using ProductManifold() constructor to generate St(d,p)^n \times Euc(p,n):
    integer numoftypes = 2;        // 2 i.e. (3D) Stiefel + Euclidean
    integer numofmani1 = somSz.n_; // num of Stiefel manifolds
    integer numofmani2 = 1;
    Stiefel mani1(somSz.p_, somSz.d_);
    mani1.ChooseParamsSet2();
    Euclidean mani2(somSz.p_, somSz.n_);
    ProductManifold ProdMani(numoftypes, &mani1, numofmani1, &mani2, numofmani2);

    // Obtain an initial iterate
    Vector startX = ProdMani.RandominManifold();

    SomUtils::readCsvInitguess("../data/rphg_test_x.csv", startX);
    startX.Print("startX");

    ProdMani.SetIsIntrApproach(false);

    // Set the domain of the problem to be the product of Stiefel manifolds
    SampleSomProblem Prob(somSz, Tijs, edges);
    Prob.SetDomain(&ProdMani);

    // Numerically check gradient consistency (optional).
    // ProdMani.CheckParams();

    Prob.SetUseGrad(true);
    Prob.SetUseHess(true);

    // output the parameters of the manifold of domain
    ROPTLIB::RTRNewton *RTRNewtonSolver = new RTRNewton(&Prob, &startX); // USE INITGUESS HERE!
    RTRNewtonSolver->Verbose = ITERRESULT;
    // RTRNewtonSolver->Max_Iteration = 500;
    // RTRNewtonSolver->Max_Inner_Iter = 500;
    // ROPTLIB::PARAMSMAP solverParams = {std::pair<std::string, double>("Max_Inner_Iter", 10)};
    // RTRNewtonSolver->SetParams(solverParams);
    RTRNewtonSolver->CheckParams();

    SomUtils::MatD xEig(SomUtils::MatD::Zero(Prob.fullSz_, 1));
    Prob.RoptToEig(startX, xEig);

    SomUtils::VecMatD R(somSz.n_, SomUtils::MatD::Zero(somSz.p_, somSz.d_));

    Prob.getRotations(xEig, R);
    SomUtils::MatD T(SomUtils::MatD::Zero(somSz.p_, somSz.n_));
    Prob.getTranslations(xEig, T);

    auto somSzNext = somSz;
    somSzNext.p_++;
    SampleSomProblem ProbNext(somSzNext, Tijs, edges);

    ROFL_VAR3(somSzNext.p_, somSzNext.d_, somSzNext.n_);
    // return;

    Stiefel mani1next(somSzNext.p_, somSz.d_);
    mani1next.ChooseParamsSet2();
    Euclidean mani2next(somSzNext.p_, somSz.n_);
    ProductManifold ProdManiNext(numoftypes, &mani1, numofmani1, &mani2, numofmani2);
    Vector Y0 = ProdManiNext.RandominManifold();
    Prob.rsomPimHessianGenproc(1e-4, R, T, Y0, false);

    Y0.Print("Y0 output of rphg small + linesearch");

    // check whether cost has actually decreased

    ROFL_VAR1(Prob.f(startX));
    ROFL_VAR1(ProbNext.f(Y0));


    delete RTRNewtonSolver;
}

int main(int argc, char **argv)
{
    int p = 3;
    int d = 3;
    int n = 10;
    SomUtils::SomSize somSz(p, d, n);
    int numEdges = 38;
    SomUtils::MatD Tijs(d, numEdges);
    Eigen::MatrixXi edges(numEdges, 2);

    SomUtils::readCsvTijs("../data/rphg_test_tijs.csv", Tijs, d, numEdges);
    SomUtils::readCsvEdges("../data/rphg_test_edges.csv", edges);

    ROFL_VAR2(Tijs.transpose(), edges); // transpose is for visualization puropses
    ROFL_VAR2(Tijs.rows(), Tijs.cols());

    testSomSample(somSz, Tijs, edges);
    return 0;
}
