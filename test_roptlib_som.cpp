
#include <iostream>
#include <fstream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <filesystem>
#include <boost/lexical_cast.hpp>

#include "thirdparty/roptlib/problems_SampleSom.h"

#include "thirdparty/roptlib/solvers_RTRNewton.h"

using namespace ROPTLIB;

void deserializeRow(const std::string &row, double &pt)
{
    std::stringstream ss(row);
    std::vector<std::string> ptCoordStrings;

    // ROFL_VAR1("\n");
    for (std::string strI; ss >> strI;)
    {
        ptCoordStrings.push_back(strI);

        if (ss.peek() == ',')
            ss.ignore();

        strI.erase(std::remove(strI.begin(), strI.end(), ','), strI.end()); // this should no nothing
        // ROFL_VAR1(strI);

        double ptCoord = std::stod(strI);

        // ROFL_VAR1(ptCoord);
        pt = ptCoord;
    }
    // ROFL_VAR1(pt.transpose());
}

void readCsvInitguess(std::string fname, Vector &csvVec)
{
    std::fstream fout;
    fout.open(fname, std::ios::in);

    if (!fout.is_open())
    {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    std::string line;
    // getline(fout, csvVec.header, '\n');
    // ROFL_VAR1(csvVec.header);
    csvVec.ObtainWriteEntireData();

    /////
    // Vector rgTiVec(sz_.p_, sz_.n_);
    realdp *GroptlibWriteArray = csvVec.ObtainWriteEntireData();

    /////

    int j = 0;
    while (getline(fout, line, '\n'))
    {
        // ROFL_VAR1(line);

        // add all the column data
        // of a row to a vector

        deserializeRow(line, GroptlibWriteArray[j]);
        j++;
        // ROFL_VAR1(line);

        // csvVec.pts.push_back(pt);
    }

    // rgTiVec.CopyTo(result->GetElement(gElemIdx));
    // csvVec.Print("csv read result");

    fout.close();
}

void testSomSample(SomSize somSz, MatD &Tijs, Eigen::MatrixXi &edges)
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

    if (false)
    {                                    // intr vs extr scope
        Vector x = startX.GetElement(0); // TODO: change later
        Vector egf(6, 1, 1);             // TODO: change later
        // egf.Initialization(6,1,1);
        Stiefel maniSt(somSz.p_, somSz.d_);
        maniSt.ChooseParamsSet1();
        Vector result(6, 1, 1);

        result = maniSt.ObtainIntr(x, egf, &result);

        // Vector tmp(somSz.d_, somSz.d_); tmp.AlphaABaddBetaThis(1, egf, GLOBAL::T, x, GLOBAL::N, 0);
        // result.AlphaABaddBetaThis(-1, x, GLOBAL::N, tmp, GLOBAL::N, 1); /*egf - x * (egf.GetTranspose() * x)*/
        result.Print("Minimum notation debugging 1");
    }

    readCsvInitguess("../data/X_initguess.csv", startX);
    startX.Print("startX");

    Vector etaStartX = ProdMani.RandominManifold();
    readCsvInitguess("../data/etaX_initguess.csv", etaStartX);
    etaStartX.Print("etaStartX");

    // ProdMani.SetIsIntrApproach(false);

    // Set the domain of the problem to be the product of Stiefel manifolds
    SampleSomProblem Prob(somSz, Tijs, edges);
    Prob.SetDomain(&ProdMani);

    // Numerically check gradient consistency (optional).
    // ProdMani.CheckParams();

    Prob.SetUseGrad(true);
    Prob.SetUseHess(true);
    // Prob.SetNumGradHess(true);
    // Prob.CheckGradHessian(startX);

    // cost f
    // realdp costStart = Prob.f(startX);
    // ROFL_VAR1(costStart);
    // grad
    // Vector gradStart = ProdMani.RandominManifold();
    // gradStart = Prob.RieGrad(startX, &gradStart);
    // gradStart.Print("grad Start");
    // hess
    Vector hessStart = ProdMani.RandominManifold();
    hessStart = Prob.RieHessianEta(startX, etaStartX, &hessStart);
    hessStart.Print("hess Start");

    // output the parameters of the manifold of domain
    ROPTLIB::RTRNewton *RTRNewtonSolver = new RTRNewton(&Prob, &startX); // USE INITGUESS HERE!
    RTRNewtonSolver->Verbose = ITERRESULT;
    // RTRNewtonSolver->Max_Iteration = 500;
    // RTRNewtonSolver->Max_Inner_Iter = 500;
    // ROPTLIB::PARAMSMAP solverParams = {std::pair<std::string, double>("Max_Inner_Iter", 10)};
    // RTRNewtonSolver->SetParams(solverParams);
    RTRNewtonSolver->CheckParams();

    // % Solve.
    // [x, xcost, info, options] = trustregions(problem);
    RTRNewtonSolver->Run();
    Prob.CheckGradHessian(RTRNewtonSolver->GetXopt());

    std::cout << "Prob.GetUseGrad() " << Prob.GetUseGrad() << std::endl;
    std::cout << "Prob.GetUseHess() " << Prob.GetUseHess() << std::endl;
    // std::cout << "Prob.GetNumGradHess() " << Prob.GetNumGradHess() << std::endl;

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
    int numEdges = 18;
    MatD Tijs(d, numEdges);
    Tijs << -7.02929319733264, -4.70228201833979, 4.34434211283999, 7.60845213036123, -7.49784261659683, -4.70228201833979, 4.52331206558989, 4.70228201833979, 7.13990271109704, -7.31887266384694, -4.70228201833979, -4.70228201833978, 4.70228201833978, 7.31887266384693, -7.13990271109703, -4.52331206558989, 4.70228201833979, 7.49784261659683,
        3.06146745892072, 0, 3.06146745892072, -0.946045472532388, -0.946045472532388, 8.88178419700125e-16, -0.946045472532388, 0, -2.47677920199275, -2.47677920199275, 0, 2.47677920199275, -8.88178419700125e-16, 2.47677920199275, 2.47677920199275, 0.946045472532387, 0, 0.946045472532388,
        5.71604418959064, 14.4721359549996, 13.9794739404545, 5.71604418959065, 5.37562311002300, 14.4721359549996, 14.5302868176819, 14.4721359549996, 5.37562311002299, 5.92643598725038, 14.4721359549996, 13.9794739404545, 14.4721359549996, 5.92643598725039, 5.37562311002299, 14.5302868176819, 14.4721359549996, 5.37562311002299;
    Eigen::MatrixXi edges(numEdges, 2);
    edges << 2, 1,
        3, 1,
        4, 1,
        1, 2,
        3, 2,
        4, 2,
        5, 2,
        1, 3,
        2, 3,
        4, 3,
        5, 3,
        1, 4,
        2, 4,
        3, 4,
        5, 4,
        2, 5,
        3, 5,
        4, 5;

    testSomSample(somSz, Tijs, edges);
    return 0;
}
