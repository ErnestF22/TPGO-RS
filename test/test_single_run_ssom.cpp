#include <iostream>
#include <fstream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <filesystem>
#include <set>
#include <boost/lexical_cast.hpp>

#include "thirdparty/roptlib/problems_Ssom.h"

#include <rofl/common/param_map.h>
#include <rofl/common/profiler.h>

namespace fs = std::filesystem;

void EigToRopt(const SomUtils::MatD &xEig, const ROPTLIB::SsomProblem &Prob, ROPTLIB::Vector *result)
{
    int rotSz = Prob.getRotSz();
    int translSz = Prob.getTranslSz();

    auto rgR = SomUtils::VecMatD(Prob.sz_.n_, SomUtils::MatD::Zero(Prob.sz_.p_, Prob.sz_.d_));
    // ROFL_VAR3(xEig.rows(), R.size(), T.size());
    // ROFL_VAR3(R.size(), P.size(), rgR.size());
    Prob.getRotations(xEig, rgR);
    SomUtils::MatD rgT(SomUtils::MatD::Zero(Prob.sz_.p_, Prob.sz_.n_));
    Prob.getTranslations(xEig, rgT);
    SomUtils::MatD rgLambdas(SomUtils::MatD::Zero(Prob.numEdges_, 1));
    Prob.getScales(xEig, rgLambdas);

    auto sz = Prob.sz_;

    int gElemIdx = 0;
    // fill result with computed gradient values: R
    for (int i = 0; i < sz.n_; ++i)
    {
        // ROFL_VAR1(gElemIdx);
        // ROFL_VAR2("\n", rgR[gElemIdx]);
        // result->GetElement(gElemIdx).SetToIdentity(); // Ri
        // result->GetElement(gElemIdx).Print("Ri before assignment");

        ROPTLIB::Vector rgRiVec(sz.p_, sz.d_);
        // rgRiVec.Initialize();
        realdp *GroptlibWriteArray = rgRiVec.ObtainWriteEntireData();
        for (int j = 0; j < rotSz; ++j)
        {
            // ROFL_VAR2(i, j);
            // rgRiVec.Print("rgRiVec before assignment");

            // ROFL_VAR1(rgRiVec.GetElement(j, 0));

            GroptlibWriteArray[j] = rgR[i].reshaped(sz.d_ * sz.p_, 1)(j);

            // ROFL_VAR1("");
            // rgRiVec.Print("rgRiVec after assignment");
        }
        rgRiVec.CopyTo(result->GetElement(gElemIdx));
        // result->GetElement(gElemIdx).Print("Riem. grad Ri after assignment");
        gElemIdx++;
    }

    // fill result with computed gradient values: T

    ROPTLIB::Vector rgTiVec(sz.p_, sz.n_);
    realdp *GroptlibWriteArray = rgTiVec.ObtainWriteEntireData();
    for (int j = 0; j < sz.p_ * sz.n_; ++j)
    {
        // rgTiVec.Print("rgTiVec before assignment");

        // ROFL_VAR1(rgRiVec.GetElement(j, 0));

        GroptlibWriteArray[j] = rgT.reshaped(sz.n_ * sz.p_, 1)(j);

        // ROFL_VAR1("");
        // rgTiVec.Print("rgTiVec after assignment");
    }
    rgTiVec.CopyTo(result->GetElement(gElemIdx));
    gElemIdx++;
    // result->GetElement(gElemIdx).Print("grad Ti after assignment");

    // ROFL_VAR2("\n", rgT);

    // fill result with computed gradient values: Lambdas

    ROPTLIB::Vector rgLambdasIvec(Prob.numEdges_, 1);
    realdp *GroptlibWriteArray2 = rgLambdasIvec.ObtainWriteEntireData();
    for (int j = 0; j < Prob.numEdges_; ++j)
    {
        // rhTiVec.Print("rhTiVec before assignment");

        // ROFL_VAR1(rhRiVec.GetElement(j, 0));

        GroptlibWriteArray2[j] = rgLambdas(j); // TODO: reshaped() call can probably be removed

        // ROFL_VAR1("");
        // rhTiVec.Print("rhTiVec after assignment");
    }
    rgLambdasIvec.CopyTo(result->GetElement(gElemIdx));
}

int main(int argc, char **argv)
{
    std::string filenameCfg;
    std::string folderIn;
    std::set<fs::path> sortedByName;

    int d;
    // int numTestsPerInstance;
    bool readStartingPtFromFile;
    int srcNodeIdx;
    double rho;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>(
        "in", folderIn,
        std::string("../matlab/data/ssom_testdata_noisy/harder/tdata_n5_mindeg2_sigma00/"));

    params.getParam<int>("d", d, 3);
    params.getParam<double>("rho", rho, 1000.0);
    // params.getParam<int>("numTestsPerInstance", numTestsPerInstance, 30);
    params.getParam<bool>("readStartingPtFromFile", readStartingPtFromFile, true);
    params.getParam<int>("srcNodeIdx", srcNodeIdx, 0);

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    std::cout << "-------\n"
              << std::endl;

    /*************************End of ROFL params reading**************************/

    int n;
    if (!SomUtils::readSingleIntCsv(folderIn + "n.csv", n))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }

    int numEdges;
    if (!SomUtils::readSingleIntCsv(folderIn + "e.csv", numEdges))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }

    ROFL_VAR2(n, numEdges); // n = 10, numEdges = 48

    SomUtils::SomSize somSzD(d, d, n);
    SomUtils::MatD Tijs(d, numEdges);
    Eigen::MatrixXi edges(numEdges, 2);

    if (!SomUtils::readMatlabCsvTijs(folderIn + "tijs.csv", Tijs, d, numEdges))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }
    if (!SomUtils::readMatlabCsvEdges(folderIn + "edges.csv", edges))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }

    ROFL_VAR1(Tijs)
    ROFL_VAR1(edges)

    // problem d x d x n
    integer numoftypes = 3; // 2 i.e. (3D) Stiefel + Euclidean
    integer numofmani1 = n; // num of Stiefel manifolds
    integer numofmani2 = 1;
    integer numofmani3 = 1;

    ROPTLIB::Stiefel mani1(d, d);
    mani1.ChooseParamsSet2();
    ROPTLIB::Euclidean mani2(d, n);
    ROPTLIB::Euclidean mani3(numEdges);
    ROPTLIB::ProductManifold ProdManiSsom(numoftypes,
                                          &mani1, numofmani1, &mani2, numofmani2, &mani3, numofmani3);
    ROPTLIB::SsomProblem Prob(somSzD, Tijs, edges);
    Prob.setRho(rho); // default 1000.0

    // Read GT from csv
    ROPTLIB::Vector xGt = ProdManiSsom.RandominManifold();
    if (!SomUtils::readCsvInitguess(folderIn + "Xgt.csv", xGt))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }

    xGt.Print("xGt");
    // ROPT to Eig (GT)
    SomUtils::MatD XgtVecEig(SomUtils::MatD::Zero(d * d * n + d * n + numEdges, 1));
    SomUtils::VecMatD RgtEig(n, SomUtils::MatD::Zero(d, d));
    SomUtils::MatD TgtEig(SomUtils::MatD::Zero(d, n));
    SomUtils::MatD LambdasGtEig(SomUtils::MatD::Zero(numEdges, 1));
    Prob.RoptToEig(xGt, XgtVecEig);
    Prob.getRotations(XgtVecEig, RgtEig);
    Prob.getTranslations(XgtVecEig, TgtEig);
    Prob.getScales(XgtVecEig, LambdasGtEig);
    Prob.setGt(RgtEig, TgtEig, LambdasGtEig);

    // Set the domain of the problem to be the product of Stiefel manifolds
    Prob.SetDomain(&ProdManiSsom);

    // Set Prob params
    Prob.SetUseGrad(true);
    Prob.SetUseHess(true);

    ROFL_VAR1("Printing R, T, Lambdas gt")
    for (auto &m : RgtEig)
        ROFL_VAR1(m)
    ROFL_VAR1(TgtEig)
    ROFL_VAR1(LambdasGtEig)

    ROFL_VAR1(Prob.costEigen(RgtEig, TgtEig, LambdasGtEig));

    // // problem nrs x d x n
    // ROPTLIB::Stiefel mani1nrs(nrs, d);
    // mani1.ChooseParamsSet2();
    // ROPTLIB::Euclidean mani2nrs(nrs, n);
    // ROPTLIB::ProductManifold ProdManiSsomNrs(numoftypes, &mani1nrs, numofmani1, &mani2nrs, numofmani2);
    // // Read from csv
    // ROPTLIB::Vector xManoptOut = ProdManiSsomNrs.RandominManifold();
    // SomUtils::readCsvInitguess("../data/recov_x_manopt_out.csv", xManoptOut);
    // xManoptOut.Print("xManoptOut");

    // Generate startX (random)
    ROPTLIB::Vector startX = ProdManiSsom.RandominManifold();

    if (readStartingPtFromFile)
        if (!SomUtils::readCsvInitguess(folderIn + "ssom_x_start.csv", startX))
        {
            // matlab/data/ssom_testdata_noisy/harder/tdata_n5_mindeg2_sigma00/ssom_x_start.csv
            ROFL_ERR("Error opening file")
            ROFL_VAR1(folderIn + "ssom_x_start.csv")
            ROFL_ASSERT(0)
        }
    startX.Print("startX");

    SomUtils::MatD startXeig = SomUtils::MatD::Zero(d * d * n + d * n + numEdges, 1);
    Prob.RoptToEig(startX, startXeig);
    SomUtils::MatD scalesInitguess = SomUtils::MatD::Ones(numEdges, 1);
    if (!Prob.reluScaleCompensation_)
    {
        scalesInitguess *= 10.0;
    }
    startXeig.block(d * d * n + d * n, 0, numEdges, 1) = scalesInitguess;

    ROPTLIB::Vector startX2 = ProdManiSsom.RandominManifold();
    EigToRopt(startXeig, Prob, &startX2);

    startX2.Print("startX2");

    SomUtils::MatD XstartVecEig(SomUtils::MatD::Zero(d * d * n + d * n + numEdges, 1));
    SomUtils::VecMatD RstartEig(n, SomUtils::MatD::Zero(d, d));
    SomUtils::MatD TstartEig(SomUtils::MatD::Zero(d, n));
    SomUtils::MatD LambdasStartEig(SomUtils::MatD::Zero(numEdges, 1));
    Prob.RoptToEig(startX2, XstartVecEig);
    Prob.getRotations(XstartVecEig, RstartEig);
    Prob.getTranslations(XstartVecEig, TstartEig);
    Prob.getScales(XstartVecEig, LambdasStartEig);

    ROFL_VAR1(Prob.costEigen(RstartEig, TstartEig, LambdasStartEig));

    ROPTLIB::Vector riegradIG = ProdManiSsom.RandominManifold();
    Prob.RieGrad(startX2, &riegradIG);
    riegradIG.Print("riegradIG");

    // ROPTLIB::Vector startU = ProdManiSsom.RandominManifold();
    // if (!SomUtils::readCsvInitguess(folderIn + "ssom_u_start.csv", startU))
    // {
    //     // matlab/data/ssom_testdata_noisy/harder/tdata_n5_mindeg2_sigma00/ssom_x_start.csv
    //     ROFL_ERR("Error opening file")
    //     ROFL_VAR1(folderIn + "ssom_u_start.csv")
    //     ROFL_ASSERT(0)
    // }
    // startU.Print("startU");

    // ROPTLIB::Vector rgradStart = ProdManiSsom.RandominManifold();
    // Prob.RieGrad(startX, &rgradStart);
    // rgradStart.Print("rgradStart");

    // SomUtils::VecMatD xR(n, SomUtils::MatD::Identity(d, d));
    // SomUtils::MatD xT(SomUtils::MatD::Zero(d, n));
    // SomUtils::MatD xLambdas(SomUtils::MatD::Zero(numEdges, 1));
    // SomUtils::MatD startXeig(SomUtils::MatD::Zero(n * d * d + d * n + numEdges, 1));
    // Prob.RoptToEig(startX, startXeig);
    // Prob.getRotations(startXeig, xR);
    // Prob.getTranslations(startXeig, xT);
    // Prob.getScales(startXeig, xLambdas);
    // SomUtils::VecMatD uR(n, SomUtils::MatD::Identity(d, d));
    // SomUtils::MatD uT(SomUtils::MatD::Zero(d, n));
    // SomUtils::MatD uLambdas(SomUtils::MatD::Zero(numEdges, 1));
    // SomUtils::MatD startUeig(SomUtils::MatD::Zero(n * d * d + d * n + numEdges, 1));
    // Prob.RoptToEig(startU, startUeig);
    // Prob.getRotations(startUeig, uR);
    // Prob.getTranslations(startUeig, uT);
    // Prob.getScales(startUeig, uLambdas);
    // SomUtils::VecMatD hessRout(n, SomUtils::MatD::Identity(d, d));
    // SomUtils::MatD hessTout(SomUtils::MatD::Zero(d, n));
    // SomUtils::MatD hessLambdasOut(SomUtils::MatD::Zero(numEdges, 1));
    // Prob.hessGenprocEigen(xR, uR, xT, uT, xLambdas, uLambdas, hessRout, hessTout, hessLambdasOut);
    // for (auto &m : hessRout)
    //     ROFL_VAR1(m)
    // ROFL_VAR1(hessTout)
    // ROFL_VAR1(hessLambdasOut)

    // RUN SSOM RS
    int srcNodeId = 0;
    SomUtils::VecMatD Rout(n, SomUtils::MatD::Identity(d, d));
    SomUtils::MatD Tout(SomUtils::MatD::Zero(d, n));
    SomUtils::MatD lambdasOut(SomUtils::MatD::Zero(numEdges, 1));
    int lastStaircaseStep;
    bool rsSuccess = false, rotDetsOk = false, lambdasAcceptable = false;
    double exectime = 0;
    {
        /* Setting up Prob using setters */
        Prob.setUsePIM(true);           // same as default
        Prob.setPimMaxIterations(5000); // same as default

        rofl::ScopedTimer timer("ssomRS");

        double costOut = ROPTLIB::runSsom(Prob, startX2, srcNodeId,
                                          Rout, Tout, lambdasOut,
                                          lastStaircaseStep,
                                          rsSuccess, rotDetsOk, lambdasAcceptable); // note: startX is needed (even if random) in ROPTLIB;
        // ROPTLIB namespace is used even if runRsomRS() is not in SsomProblem class, nor in "original" ROPTLIB
        ROFL_VAR1(costOut)
        exectime = timer.elapsedTimeMs();
    }

    std::vector<double> rotErrs(numEdges, 1e+6), translErrs(numEdges, 1e+6), scaleErrs(numEdges, 1e+6);
    SomUtils::computeErrorsSingleSsom(edges,
                                      Rout, Tout, lambdasOut,
                                      RgtEig, TgtEig, LambdasGtEig,
                                      rotErrs, translErrs, scaleErrs);

    ROFL_VAR1(exectime)

    return 0;
}
