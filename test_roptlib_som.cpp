
#include <iostream>
#include <fstream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <filesystem>
#include <boost/lexical_cast.hpp>

#include "thirdparty/roptlib/manifolds_Euclidean.h"
#include "thirdparty/roptlib/manifolds_MultiManifolds.h"

#include "thirdparty/roptlib/problems_SampleSom.h"

#include "thirdparty/roptlib/solvers_RTRNewton.h"

using namespace ROPTLIB;

void testSomSample(SomSize somSz, Eigen::MatrixXd &Tijs, Eigen::MatrixXi &edges)
{
    ROFL_VAR3(somSz.p_, somSz.d_, somSz.n_);

    // manifold
    // An example of using ProductManifold() constructor to generate St(d,p)^n \times Euc(p,n):
    integer numoftypes = 2; // 2 i.e. (3D) Stiefel + Euclidean
    integer numofmani1 = somSz.n_; // num of Stiefel manifolds
    integer numofmani2 = 1;
    Stiefel mani1(somSz.p_, somSz.d_);
    Euclidean mani2(somSz.p_, somSz.n_);
    ProductManifold ProdMani(numoftypes, &mani1, numofmani1, &mani2, numofmani2);
    

    // Obtain an initial iterate
    Vector startX = ProdMani.RandominManifold();
    startX.ObtainWriteEntireData();
    realdp *startXWriteArray = startX.ObtainWriteEntireData(); //!! assignment in col-major order
    // startXWriteArray[0] = ...;
    startX.Print("startX");

    // Set the domain of the problem to be the product of Stiefel manifolds
    SampleSomProblem Prob(somSz, Tijs, edges);
    Prob.SetDomain(&ProdMani);

    // Numerically check gradient consistency (optional).
    // ProdMani.CheckParams();
    
    Prob.SetUseGrad(false);
    // Prob.SetUseHess(false);
    Prob.SetNumGradHess(true);
    // Prob.CheckGradHessian(ProdX);

    // output the parameters of the manifold of domain
    ROPTLIB::RTRNewton *RTRNewtonSolver = new RTRNewton(&Prob, &startX); // USE INITGUESS HERE!
    RTRNewtonSolver->Verbose = ITERRESULT;
    // RTRNewtonSolver->Max_Iteration = 500;
    // RTRNewtonSolver->Max_Inner_Iter = 500;
    RTRNewtonSolver->CheckParams();
    // auto solverParams = {std::pair<std::string, double>("Max_Inner_Iter", 10)};
    // RTRNewtonSolver->SetParams(solverParams);

    // % Solve.
    // [x, xcost, info, options] = trustregions(problem);

    RTRNewtonSolver->Run();
    // Prob.CheckGradHessian(RTRNewtonSolver->GetXopt());

    std::cout << "Prob.GetUseGrad() " << Prob.GetUseGrad() << std::endl;
    std::cout << "Prob.GetUseHess() " << Prob.GetUseHess() << std::endl;
    std::cout << "Prob.GetNumGradHess() " << Prob.GetNumGradHess() << std::endl;

    // Outputs
    RTRNewtonSolver->GetXopt().Print("Xopt");
    std::cout << "xcost " << RTRNewtonSolver->Getfinalfun() << std::endl; // x cost

    delete RTRNewtonSolver;
}

int main(int argc, char **argv)
{
    int p = 4;
    int d = 3;
    int n = 5;
    SomSize somSz(p, d, n);
    int numEdges = 14;
    Eigen::MatrixXd Tijs(d, numEdges);
    Tijs << -7.02929319733264, -4.70228201833979, 7.60845213036123, -7.49784261659683, -4.70228201833979, 4.70228201833979, 7.13990271109704, -7.31887266384694, -4.70228201833979, 4.70228201833978, 7.31887266384693, -7.13990271109703, 4.70228201833979, 7.49784261659683,
        3.06146745892072, 0, -0.946045472532388, -0.946045472532388, 8.88178419700125e-16, 0, -2.47677920199275, -2.47677920199275, 0, -8.88178419700125e-16, 2.47677920199275, 2.47677920199275, 0, 0.946045472532388,
        5.71604418959064, 14.4721359549996, 5.71604418959065, 5.37562311002300, 14.4721359549996, 14.4721359549996, 5.37562311002299, 5.92643598725038, 14.4721359549996, 14.4721359549996, 5.92643598725039, 5.37562311002299, 14.4721359549996, 5.37562311002299;
    Eigen::MatrixXi edges(numEdges, 2);
    edges << 2, 1,
        3, 1,
        1, 2,
        3, 2,
        4, 2,
        1, 3,
        2, 3,
        4, 3,
        5, 3,
        2, 4,
        3, 4,
        5, 4,
        3, 5,
        4, 5;

    testSomSample(somSz, Tijs, edges);
    return 0;
}
