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

int main(int argc, char **argv)
{
    std::string filenameCfg;
    std::string folderIn;
    std::set<fs::path> sortedByName;

    int d;
    // int numTestsPerInstance;
    bool readStartingPtFromFile;
    int srcNodeIdx;

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

    // RUN RSOM RS
    int srcNodeId = 0;
    SomUtils::VecMatD Rout(n, SomUtils::MatD::Identity(d, d));
    SomUtils::MatD Tout(SomUtils::MatD::Zero(d, n));
    SomUtils::MatD lambdasOut(SomUtils::MatD::Zero(numEdges, 1));
    int lastStaircaseStep;
    double exectime = 0;
    {
        rofl::ScopedTimer timer("RsomRS");
        double costOut = ROPTLIB::runSsom(Prob, startX, srcNodeId,
                                          Rout, Tout, lambdasOut,
                                          lastStaircaseStep); // note: startX is needed (even if random) in ROPTLIB;
        // ROPTLIB namespace is used even if runRsomRS() is not in SsomProblem class, nor in "original" ROPTLIB
        ROFL_VAR1(costOut)
        exectime = timer.elapsedTimeMs();
    }

    // std::vector<double> rotErrs(numEdges, 1e+6), translErrs(numEdges, 1e+6);
    // SomUtils::computeErrorsSingleRsom(edges,
    //                                   Rout, Tout,
    //                                   RgtEig, TgtEig,
    //                                   rotErrs, translErrs);
    // ROFL_VAR1(exectime)

    return 0;
}
