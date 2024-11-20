
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

    // SomUtils::readCsvInitguess("../data/X_initguess.csv", startX);
    // startX.Print("startX");

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

    // RS step !!!!

    int nrsNext = somSz.p_ + 1;
    SomUtils::SomSize szNext(nrsNext, somSz.d_, somSz.n_);

    Vector xOptPrevStep = RTRNewtonSolver->GetXopt();
    int fullSz = somSz.p_ * somSz.d_ * somSz.n_ + somSz.p_ * somSz.n_;
    int fullSzNext = szNext.p_ * szNext.d_ * szNext.n_ + szNext.p_ * szNext.n_;

    Stiefel mani1Next(szNext.p_, szNext.d_);
    mani1Next.ChooseParamsSet2();

    Euclidean mani2Next(szNext.p_, szNext.n_);
    ProductManifold ProdManiNext(numoftypes, &mani1Next, numofmani1, &mani2Next, numofmani2);

    Vector xOptNext = ProdManiNext.RandominManifold(); //!! in other cases xOptNext would have been a pointer

    Vector Y0 = ProdManiNext.RandominManifold();

    SomUtils::MatD xOptPrevStepEig(SomUtils::MatD::Zero(fullSz, 1));

    Prob.RoptToEig(xOptPrevStep, xOptPrevStepEig);

    SomUtils::VecMatD R(somSz.n_, SomUtils::MatD::Zero(somSz.p_, somSz.d_));
    Prob.getRotations(xOptPrevStepEig, R);

    SomUtils::MatD T(SomUtils::MatD::Zero(somSz.p_, somSz.n_));
    Prob.getTranslations(xOptPrevStepEig, T);

    // 1) cat zero rows on R, T
    SomUtils::VecMatD Rnext(szNext.n_, SomUtils::MatD::Zero(szNext.p_, szNext.d_));
    Prob.catZeroRow3dArray(R, Rnext); // does the same job as ProbNext.catZeroRow3dArray()

    SomUtils::MatD Tnext(SomUtils::MatD::Zero(szNext.p_, szNext.n_));
    Prob.catZeroRow(T, Tnext); // does the same job as ProbNext.catZeroRow()

    SampleSomProblem ProbNext(szNext, Tijs, edges);
    ProbNext.SetDomain(&ProdMani);

    // 2) Eig to ROPT
    { // EigToRopt scope for xOptNext

        int rotSz = ProbNext.getRotSz();
        int translSz = ProbNext.getTranslSz();

        int gElemIdx = 0;
        // fill result with computed gradient values : R
        for (int i = 0; i < szNext.n_; ++i)
        {
            // ROFL_VAR1(gElemIdx);
            // ROFL_VAR2("\n", rgR[gElemIdx]);
            // result->GetElement(gElemIdx).SetToIdentity(); // Ri
            // result->GetElement(gElemIdx).Print("Ri before assignment");

            Vector RnextROPT(szNext.p_, szNext.d_);
            // RnextROPT.Initialize();
            realdp *GroptlibWriteArray = RnextROPT.ObtainWriteEntireData();
            for (int j = 0; j < rotSz; ++j)
            {
                // ROFL_VAR2(i, j);
                // RnextROPT.Print("RnextROPT before assignment");

                // ROFL_VAR1(RnextROPT.GetElement(j, 0));

                GroptlibWriteArray[j] = Rnext[i].reshaped(szNext.d_ * szNext.p_, 1)(j);

                // ROFL_VAR1("");
                // RnextROPT.Print("RnextROPT after assignment");
            }
            RnextROPT.CopyTo(xOptNext.GetElement(gElemIdx));
            // result->GetElement(gElemIdx).Print("Riem. grad Ri after assignment");
            gElemIdx++;
        }

        // fill result with computed gradient values : T

        Vector TnextROPT(szNext.p_, szNext.n_);
        realdp *GroptlibWriteArray = TnextROPT.ObtainWriteEntireData();
        for (int j = 0; j < szNext.p_ * szNext.n_; ++j)
        {
            // TnextROPT.Print("TnextROPT before assignment");

            // ROFL_VAR1(RnextROPT.GetElement(j, 0));

            GroptlibWriteArray[j] = Tnext.reshaped(szNext.n_ * szNext.p_, 1)(j);

            // ROFL_VAR1("");
            // TnextROPT.Print("TnextROPT after assignment");
        }
        TnextROPT.CopyTo(xOptNext.GetElement(gElemIdx));
    } // end of EigToRopt scope for xOptNext

    //

    Prob.linesearchArmijoROPTLIB(xOptNext, szNext, Y0);

    Y0.Print("Y0 before delete RTRNewtonSolver");

    delete RTRNewtonSolver;
}

int main(int argc, char **argv)
{
    int p = 4;
    int d = 3;
    int n = 5;
    SomUtils::SomSize somSz(p, d, n);
    int numEdges = 18;
    SomUtils::MatD Tijs = SomUtils::MatD::Random(d, numEdges);
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
