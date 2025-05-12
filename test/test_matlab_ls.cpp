#include <iostream>
#include <fstream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <filesystem>
#include <boost/lexical_cast.hpp>

#include "thirdparty/qr_unique_sizeless/main.h"

#include "thirdparty/roptlib/problems_SampleSom.h"

int main(int argc, char **argv)
{
    int p = 4;
    int d = 3;
    int n;
    SomUtils::readSingleIntCsv("../matlab/data/cpp_testdata/tdata_n10_mindeg3/n.csv", n);

    int numEdges;

    SomUtils::readSingleIntCsv("../matlab/data/cpp_testdata/tdata_n10_mindeg3/e.csv", numEdges);

    ROFL_VAR2(n, numEdges); // n = 10, numEdges = 48

    SomUtils::SomSize somSz(p, d, n);
    SomUtils::MatD Tijs(d, numEdges);
    Eigen::MatrixXi edges(numEdges, 2);

    SomUtils::readMatlabCsvTijs("../matlab/data/cpp_testdata/tdata_n10_mindeg3/tijs.csv", Tijs, d, numEdges);
    SomUtils::readMatlabCsvEdges("../matlab/data/cpp_testdata/tdata_n10_mindeg3/edges.csv", edges);

    ROFL_VAR1(Tijs)
    ROFL_VAR1(edges)

    // problem p x d x n
    integer numoftypes = 2; // 2 i.e. (3D) Stiefel + Euclidean
    integer numofmani1 = n; // num of Stiefel manifolds
    integer numofmani2 = 1;
    ROPTLIB::Stiefel mani1(p, d);
    mani1.ChooseParamsSet2();
    ROPTLIB::Euclidean mani2(p, n);
    ROPTLIB::ProductManifold ProdMani(numoftypes, &mani1, numofmani1, &mani2, numofmani2);
    ROPTLIB::SampleSomProblem Prob(somSz, Tijs, edges);

    // Read Xnext from csv
    ROPTLIB::Vector xNext = ProdMani.RandominManifold();
    SomUtils::readCsvInitguess("../matlab/data/lsdummy_debug/matlab_Xnext_vec.csv", xNext);
    xNext.Print("xNext");
    // ROPT to Eig (Xnext)
    SomUtils::MatD xNextVecEig(SomUtils::MatD::Zero(p * d * n + p * n, 1));
    SomUtils::VecMatD RnextEig(n, SomUtils::MatD::Zero(p, d));
    SomUtils::MatD TnextEig(SomUtils::MatD::Zero(p, n));
    Prob.RoptToEig(xNext, xNextVecEig);
    ROFL_VAR2(xNextVecEig.rows(), xNextVecEig.cols())
    Prob.getRotations(xNextVecEig, RnextEig);
    Prob.getTranslations(xNextVecEig, TnextEig);
    ROFL_VAR1(Prob.f(xNext));

    // Read vpas from csv
    ROPTLIB::Vector vpas = ProdMani.RandominManifold();
    SomUtils::readCsvInitguess("../matlab/data/lsdummy_debug/matlab_vpas_vec.csv", vpas);
    vpas.Print("vpas");
    // ROPT to Eig (vpas)
    SomUtils::MatD vpasVecEig(SomUtils::MatD::Zero(p * d * n + p * n, 1));
    SomUtils::VecMatD vpasREig(n, SomUtils::MatD::Zero(p, d));
    SomUtils::MatD vpasTEig(SomUtils::MatD::Zero(p, n));
    Prob.RoptToEig(vpas, vpasVecEig);
    ROFL_VAR2(vpasVecEig.rows(), vpasVecEig.cols())
    Prob.getRotations(vpasVecEig, vpasREig);
    Prob.getTranslations(vpasVecEig, vpasTEig);

    // Read Y0 from csv
    ROPTLIB::Vector Y0mlab = ProdMani.RandominManifold();
    SomUtils::readCsvInitguess("../matlab/data/lsdummy_debug/matlab_Y0_vec.csv", xNext);
    xNext.Print("Y0mlab");
    // ROPT to Eig (Y0mlab)
    SomUtils::MatD Y0mlabVecEig(SomUtils::MatD::Zero(p * d * n + p * n, 1));
    SomUtils::VecMatD Y0mlabREig(n, SomUtils::MatD::Zero(p, d));
    SomUtils::MatD Y0mlabTEig(SomUtils::MatD::Zero(p, n));
    Prob.RoptToEig(Y0mlab, Y0mlabVecEig);
    ROFL_VAR2(Y0mlabVecEig.rows(), Y0mlabVecEig.cols())
    Prob.getRotations(Y0mlabVecEig, Y0mlabREig);
    Prob.getTranslations(Y0mlabVecEig, Y0mlabTEig);

    double costBeforeLS;
    SomUtils::readSingleDoubleCsv("../matlab/data/lsdummy_debug/matlab_cost_before_ls.csv", costBeforeLS);

    double costAfterLS;
    SomUtils::readSingleDoubleCsv("../matlab/data/lsdummy_debug/matlab_cost_after_ls.csv", costAfterLS);

    // Run linesearch dummy
    SomUtils::VecMatD Y0R(n, SomUtils::MatD::Zero(p, d));
    SomUtils::MatD Y0T(SomUtils::MatD::Zero(p, n));
    bool isQRretr = true;
    Prob.linesearchDummy(costBeforeLS, RnextEig, TnextEig, vpasREig, vpasTEig, Y0R, Y0T, isQRretr);

    ROFL_VAR2(costBeforeLS, Prob.costEigen(Y0R, Y0T));

    // Debug also RSOM Pim Hessian Genproc -> seems OK -> commented code below
    // SomUtils::VecMatD Reig(n, SomUtils::MatD::Zero(p, d));
    // for (int i = 0; i < n; ++i)
    //     Reig[i] = RnextEig[i].block(0, 0, d, d);

    // SomUtils::MatD Teig(SomUtils::MatD::Zero(d, n));
    // Teig = TnextEig.block(0, 0, d, n);
    // ROPTLIB::Vector Y0cpp = ProdMani.RandominManifold();
    // Prob.rsomPimHessianGenproc(1e-6, Reig, Teig, Y0cpp);

    /**
     * Try optimization starting from 2 points
     */
    SomUtils::SomSize somSzNext(4, d, n);

    ROPTLIB::SampleSomProblem ProbNext(somSzNext, Tijs, edges);
    ROPTLIB::Stiefel mani1next(somSzNext.p_, somSzNext.d_);
    mani1next.ChooseParamsSet2();
    ROPTLIB::Euclidean mani2next(somSzNext.p_, somSzNext.n_);
    ROPTLIB::ProductManifold ProdManiNext(numoftypes, &mani1next, numofmani1, &mani2next, numofmani2);
    ROPTLIB::Vector Y0 = ProdManiNext.RandominManifold();
    ProbNext.SetDomain(&ProdManiNext);
    ProbNext.SetUseGrad(true);
    ProbNext.SetUseHess(true);


    // 2) Eig to ROPT
    { // EigToRopt scope for Y0

        int rotSz = ProbNext.getRotSz();
        int translSz = ProbNext.getTranslSz();

        int gElemIdx = 0;
        // fill result with computed gradient values : R
        for (int i = 0; i < somSzNext.n_; ++i)
        {
            // ROFL_VAR1(gElemIdx);
            // ROFL_VAR2("\n", rgR[gElemIdx]);
            // result->GetElement(gElemIdx).SetToIdentity(); // Ri
            // result->GetElement(gElemIdx).Print("Ri before assignment");

            ROPTLIB::Vector Y0RROPT(somSzNext.p_, somSzNext.d_);
            // Y0RROPT.Initialize();
            realdp *GroptlibWriteArray = Y0RROPT.ObtainWriteEntireData();
            for (int j = 0; j < rotSz; ++j)
            {
                // ROFL_VAR2(i, j);
                // Y0RROPT.Print("Y0RROPT before assignment");

                // ROFL_VAR1(Y0RROPT.GetElement(j, 0));

                GroptlibWriteArray[j] = Y0R[i].reshaped(somSzNext.d_ * somSzNext.p_, 1)(j);

                // ROFL_VAR1("");
                // Y0RROPT.Print("Y0RROPT after assignment");
            }
            Y0RROPT.CopyTo(Y0.GetElement(gElemIdx));
            // result->GetElement(gElemIdx).Print("Riem. grad Ri after assignment");
            gElemIdx++;
        }

        // fill result with computed gradient values : T

        ROPTLIB::Vector Y0TROPT(somSzNext.p_, somSzNext.n_);
        realdp *GroptlibWriteArray = Y0TROPT.ObtainWriteEntireData();
        for (int j = 0; j < somSzNext.p_ * somSzNext.n_; ++j)
        {
            // Y0TROPT.Print("Y0TROPT before assignment");

            GroptlibWriteArray[j] = Y0T.reshaped(somSzNext.n_ * somSzNext.p_, 1)(j);

            // ROFL_VAR1("");
            // Y0TROPT.Print("Y0TROPT after assignment");
        }
        Y0TROPT.CopyTo(Y0.GetElement(gElemIdx));
    } // end of EigToRopt scope for Y0


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
    // ProbNext.CheckGradHessian(XoptNext);

    double costLast = XoptNextCost;

    ROFL_VAR1(costLast)
}