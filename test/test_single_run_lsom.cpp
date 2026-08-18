
#include <iostream>
#include <fstream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <filesystem>
#include <set>
#include <boost/lexical_cast.hpp>

#include "thirdparty/roptlib/problems_Lsom.h"

#include <rofl/common/param_map.h>
#include <rofl/common/profiler.h>

namespace fs = std::filesystem;

int main(int argc, char **argv)
{
    std::string filenameCfg;
    std::string folderIn;
    std::set<fs::path> sortedByName;

    int d;
    // int numTestsPerInstance;
    bool readStartingPtFromFile;
    int srcNodeIdx;
    double muAdmm;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>(
        "in", folderIn,
        std::string("../matlab/data/ssom_testdata_noisy/easy/"));

    params.getParam<int>("d", d, 3);
    // params.getParam<int>("numTestsPerInstance", numTestsPerInstance, 30);
    params.getParam<bool>("readStartingPtFromFile", readStartingPtFromFile, true);
    params.getParam<int>("srcNodeIdx", srcNodeIdx, 0);
    params.getParam<double>("muAdmm", muAdmm, 0.1);

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

    if (!SomUtils::readMatlabCsvTijs(folderIn + "tijs_truth.csv", Tijs, d, numEdges))
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
    ROPTLIB::ProductManifold ProdManiLsom(numoftypes,
                                          &mani1, numofmani1, &mani2, numofmani2, &mani3, numofmani3);

    SomUtils::MatD TijsNois = Tijs;
    for (int e = 0; e < numEdges; ++e)
    {
        TijsNois.col(e) = Tijs.col(e) / Tijs.col(e).norm();
    }
    ROPTLIB::LsomProblem Prob(somSzD, TijsNois, edges);

    // Read GT from csv
    ROPTLIB::Vector xGt = ProdManiLsom.RandominManifold();
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
    Prob.setMuAdmm(muAdmm);
    Prob.setPerformGlobalization(true);

    // Set the domain of the problem to be the product of Stiefel manifolds
    Prob.SetDomain(&ProdManiLsom);

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
    // ROPTLIB::ProductManifold ProdManiSsomNrs(numoftypes, &mani1nrs, numofmani1, &mani2nrs, numofmani2);
    // // Read from csv
    // ROPTLIB::Vector xManoptOut = ProdManiSsomNrs.RandominManifold();
    // SomUtils::readCsvInitguess("../data/recov_x_manopt_out.csv", xManoptOut);
    // xManoptOut.Print("xManoptOut");

    ROPTLIB::Vector startX = ProdManiLsom.RandominManifold();

    if (readStartingPtFromFile)
        if (!SomUtils::readCsvInitguess(folderIn + "ssom_x_start.csv", startX))
        {
            ROFL_ERR("Error opening file")
            ROFL_ASSERT(0)
        }

    startX.Print("startX");
    SomUtils::MatD startXeig(SomUtils::MatD::Zero(d * d * n + d * n + numEdges, 1));
    SomUtils::VecMatD RstartEig(n, SomUtils::MatD::Zero(d, d));
    SomUtils::MatD TstartEig(SomUtils::MatD::Zero(d, n));
    SomUtils::MatD LambdasStartEig(SomUtils::MatD::Zero(numEdges, 1));
    Prob.RoptToEig(startX, startXeig);
    Prob.getRotations(startXeig, RstartEig);
    Prob.getTranslations(startXeig, TstartEig);
    // Prob.getScales(startXeig, LambdasStartEig);

    if (Prob.reluScaleCompensation_)
        LambdasStartEig = 5.0 * SomUtils::MatD::Ones(numEdges, 1);
    else
        LambdasStartEig = 10.0 * SomUtils::MatD::Ones(numEdges, 1);

    auto startX2 = ProdManiLsom.RandominManifold();

    int gElemIdx = 0;
    // fill result with computed gradient values: R
    for (int i = 0; i < n; ++i)
    {
        // ROFL_VAR1(gElemIdx);
        // ROFL_VAR2("\n", rhR[gElemIdx]);
        // result->GetElement(gElemIdx).SetToIdentity(); // Ri
        // result->GetElement(gElemIdx).Print("Ri before assignment");

        ROPTLIB::Vector rhRiVec(d, d);
        int rotSz = d * d;
        // rhRiVec.Initialize();
        realdp *GroptlibWriteArray = rhRiVec.ObtainWriteEntireData();
        for (int j = 0; j < rotSz; ++j)
        {
            // ROFL_VAR2(i, j);
            // rhRiVec.Print("rhRiVec before assignment");

            // ROFL_VAR1(rhRiVec.GetElement(j, 0));

            GroptlibWriteArray[j] = RstartEig[i].reshaped(d * d, 1)(j);

            // ROFL_VAR1("");
            // rhRiVec.Print("rhRiVec after assignment");
        }
        rhRiVec.CopyTo(startX2.GetElement(gElemIdx));
        // result->GetElement(gElemIdx).Print("Riem. grad Ri after assignment");
        gElemIdx++;
    }

    // fill result with computed gradient values: T

    ROPTLIB::Vector rhTiVec(d, n);
    realdp *GroptlibWriteArray = rhTiVec.ObtainWriteEntireData();
    for (int j = 0; j < d * n; ++j)
    {
        // rhTiVec.Print("rhTiVec before assignment");

        // ROFL_VAR1(rhRiVec.GetElement(j, 0));

        GroptlibWriteArray[j] = TstartEig.reshaped(n * d, 1)(j);

        // ROFL_VAR1("");
        // rhTiVec.Print("rhTiVec after assignment");
    }
    rhTiVec.CopyTo(startX2.GetElement(gElemIdx));
    gElemIdx++;
    // result->GetElement(gElemIdx).Print("RieHess T after assignment");

    // fill result with computed gradient values: Lambdas

    ROPTLIB::Vector rhLambdasIvec(numEdges, 1);
    realdp *GroptlibWriteArray2 = rhLambdasIvec.ObtainWriteEntireData();
    for (int j = 0; j < numEdges; ++j)
    {
        // rhTiVec.Print("rhTiVec before assignment");

        // ROFL_VAR1(rhRiVec.GetElement(j, 0));

        GroptlibWriteArray2[j] = LambdasStartEig.reshaped(numEdges, 1)(j); // TODO: reshaped() call can probably be removed

        // ROFL_VAR1("");
        // rhTiVec.Print("rhTiVec after assignment");
    }
    rhLambdasIvec.CopyTo(startX2.GetElement(gElemIdx));

    Prob.setZAdmm(LambdasStartEig);

    // return 0;

    // RUN LSOM
    SomUtils::VecMatD Rout(n, SomUtils::MatD::Identity(d, d));
    SomUtils::MatD Tout(SomUtils::MatD::Zero(d, n));
    SomUtils::MatD lambdasOut(SomUtils::MatD::Zero(numEdges, 1));
    int lastStaircaseStep = 3;
    bool rsSuccess = false, rotDetsOk = false, lambdasAcceptable = false, rsActuallyUseful = true;
    double exectime = 0;
    ROFL_VAR1("Before ROPTLIB::runLsom")
    ROFL_VAR1(Prob.f(startX2))
    {
        rofl::ScopedTimer timer("Lsom");
        double costOut = ROPTLIB::runLsom(Prob, startX2, srcNodeIdx,
                                          Rout, Tout, lambdasOut,
                                          lastStaircaseStep,
                                          rsSuccess, rotDetsOk, lambdasAcceptable, rsActuallyUseful); // note: startX is needed (even if random) in ROPTLIB;
        // ROPTLIB namespace is used even if runLsom() is not in SsomProblem class, nor in "original" ROPTLIB
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
