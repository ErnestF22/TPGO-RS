#include <iostream>
#include <fstream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <filesystem>
#include <set>
#include <boost/lexical_cast.hpp>

#include "thirdparty/roptlib/problems_SampleSom.h"

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
        std::string("../matlab/data/cpp_testdata_noisy/tdata_n5_mindeg3_sigma00/"));

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
    integer numoftypes = 2; // 2 i.e. (3D) Stiefel + Euclidean
    integer numofmani1 = n; // num of Stiefel manifolds
    integer numofmani2 = 1;
    ROPTLIB::Stiefel mani1(d, d);
    mani1.ChooseParamsSet2();
    ROPTLIB::Euclidean mani2(d, n);
    ROPTLIB::ProductManifold ProdMani(numoftypes, &mani1, numofmani1, &mani2, numofmani2);
    ROPTLIB::SampleSomProblem Prob(somSzD, Tijs, edges);

    // Read GT from csv
    ROPTLIB::Vector xGt = ProdMani.RandominManifold();
    if (!SomUtils::readCsvInitguess(folderIn + "Xgt.csv", xGt))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }

    xGt.Print("xGt");
    // ROPT to Eig (GT)
    SomUtils::MatD XgtVecEig(SomUtils::MatD::Zero(d * d * n + d * n, 1));
    SomUtils::VecMatD RgtEig(n, SomUtils::MatD::Zero(d, d));
    SomUtils::MatD TgtEig(SomUtils::MatD::Zero(d, n));
    Prob.RoptToEig(xGt, XgtVecEig);
    Prob.getRotations(XgtVecEig, RgtEig);
    Prob.getTranslations(XgtVecEig, TgtEig);
    Prob.setGt(RgtEig, TgtEig);

    // Set the domain of the problem to be the product of Stiefel manifolds
    Prob.SetDomain(&ProdMani);

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
    // ROPTLIB::ProductManifold ProdManiNrs(numoftypes, &mani1nrs, numofmani1, &mani2nrs, numofmani2);
    // // Read from csv
    // ROPTLIB::Vector xManoptOut = ProdManiNrs.RandominManifold();
    // SomUtils::readCsvInitguess("../data/recov_x_manopt_out.csv", xManoptOut);
    // xManoptOut.Print("xManoptOut");

    // Generate startX (random)
    ROPTLIB::Vector startX = ProdMani.RandominManifold();

    if (readStartingPtFromFile)
        if (!SomUtils::readCsvInitguess(folderIn + "startX.csv", startX))
        {
            ROFL_ERR("Error opening file")
            ROFL_ASSERT(0)
        }

    // RUN RSOM RS
    int srcNodeId = 0;
    SomUtils::VecMatD Rout(n, SomUtils::MatD::Identity(d, d));
    SomUtils::MatD Tout(SomUtils::MatD::Zero(d, n));
    int lastStaircaseStep;
    double exectime = 0;
    {
        rofl::ScopedTimer timer("RsomRS");
        double costOut = ROPTLIB::runRsomRS(Prob, startX, srcNodeId, Rout, Tout, lastStaircaseStep); // note: startX is needed (even if random) in ROPTLIB;
        // ROPTLIB namespace is used even if runRsomRS() is not in SampleSomProblem class, nor in "original" ROPTLIB
        ROFL_VAR1(costOut)
        exectime = timer.elapsedTimeMs();
    }

    std::vector<double> rotErrs(numEdges, 1e+6), translErrs(numEdges, 1e+6);
    SomUtils::computeErrorsSingleRsom(edges,
                                      Rout, Tout,
                                      RgtEig, TgtEig,
                                      rotErrs, translErrs);
    ROFL_VAR1(exectime)

    return 0;
}