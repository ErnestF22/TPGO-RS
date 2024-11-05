
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

void testSomSample(SomSize somSz, Eigen::MatrixXf& Tijs, Eigen::MatrixXi& edges)
{
    integer d = 3;
    // manifold
    // Rotations manif(d);
    Stiefel manif(d, d);
    // Obtain an initial iterate
    Vector startX = manif.RandominManifold();
    startX.ObtainWriteEntireData();
    realdp *startXWriteArray = startX.ObtainWriteEntireData(); //!! assignment in col-major order
    // startXWriteArray[0] = 0.984489045287494;
    // startXWriteArray[1] = 0.175180074041903;
    // startXWriteArray[2] = -0.00965719253155387;
    // startXWriteArray[3] = -0.168366666411450;
    // startXWriteArray[4] = 0.927853513799599;
    // startXWriteArray[5] = -0.332776986240385;
    // startXWriteArray[6] = -0.0493354370651904;
    // startXWriteArray[7] = 0.329241246790877;
    // startXWriteArray[8] = 0.942956105055360;
    startX.Print("startX");

    // Set the domain of the problem to be the product of Stiefel manifolds
    SampleSomProblem Prob(somSz, Tijs, edges);
    Prob.SetDomain(&manif);

    // % Numerically check gradient consistency (optional).
    // checkgradient(problem);
    manif.CheckParams();
    Prob.SetUseGrad(true);
    Prob.SetUseHess(false);
    // Prob.SetNumGradHess(true);
    // Prob.CheckGradHessian(ProdX);

    // output the parameters of the manifold of domain
    ROPTLIB::RTRNewton *RTRNewtonSolver = new RTRNewton(&Prob, &startX); // USE INITGUESS HERE!
    RTRNewtonSolver->Verbose = ITERRESULT;
    RTRNewtonSolver->Max_Iteration = 500;
    RTRNewtonSolver->Max_Inner_Iter = 500;
    RTRNewtonSolver->CheckParams();
    // PARAMSMAP solverParams = {std::pair<std::string, double>("Max_Inner_Iter", 10)};
    // RTRNewtonSolver->SetParams(solverParams);

    // % Solve.
    // [x, xcost, info, options] = trustregions(problem);

    RTRNewtonSolver->Run();
    Prob.CheckGradHessian(RTRNewtonSolver->GetXopt());

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
    // An example of using this constructor to generate St(2,3)^2 \times Euc(2,2,3):

    integer n = 3, p = 2, m = 2;
    integer numoftypes = 2;
    integer numofmani1 = 2; // num of Stiefel manifolds
    integer numofmani2 = 1;
    Stiefel mani1(n, p);
    Euclidean mani2(m, p, n);
    ProductManifold ProdMani(numoftypes, &mani1, numofmani1, &mani2, numofmani2);

    SomSize somSz;
    Eigen::MatrixXf Tijs;
    Eigen::MatrixXi edges;

    testSomSample(somSz, Tijs, edges);
    return 0;
}
