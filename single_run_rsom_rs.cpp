#include <iostream>
#include <fstream>
#include <set>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <filesystem>
#include <boost/lexical_cast.hpp>

#include <rofl/common/param_map.h>
#include <rofl/common/profiler.h>

#include "thirdparty/roptlib/problems_SampleSom.h"

int main(int argc, char **argv)
{
    std::string filenameCfg;

    int d, n, numEdges;
    int numTestsPerInstance;
    bool readStartingPtFromFile;
    int srcNodeIdx;
    std::string basepath;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)

    params.getParam<int>("d", d, 3);
    params.getParam<int>("n", n, 5);
    params.getParam<int>("numEdges", numEdges, 18);
    // params.getParam<int>("numTestsPerInstance", numTestsPerInstance, 30);
    params.getParam<std::string>("basepath", basepath, "../results/n5_mindeg3_noise00_20250204_1223_43/n5_mindeg3_noise00_j0_");

    params.getParam<bool>("readStartingPtFromFile", readStartingPtFromFile, false);
    // params.getParam<std::string>("resultsBasePath", resultsBasePath, "../results/");
    params.getParam<int>("srcNodeIdx", srcNodeIdx, 0);

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    std::cout << "-------\n"
              << std::endl;

    /*************************End of ROFL params reading**************************/


    // "som/results/n5_mindeg3_noise00_20250204_1223_43/n5_mindeg3_noise00_j0_tijs.txt"


    ROFL_VAR2(n, numEdges); // n = 5, numEdges = 18

    SomUtils::SomSize somSzD(d, d, n);
    SomUtils::MatD Tijs(d, numEdges);
    Eigen::MatrixXi edges(numEdges, 2);

    if (!SomUtils::readMatlabCsvTijs(basepath + "tijs.txt", Tijs, d, numEdges))
    {
        ROFL_ERR("Error opening file " + basepath + "tijs.txt")
        ROFL_ASSERT(0)
    }
    if (!SomUtils::readMatlabCsvEdges(basepath + "edges.txt", edges))
    {
        ROFL_ERR("Error opening file " + basepath + "edges.txt")
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
    if (!SomUtils::readCsvInitguess("../matlab/data/cpp_testdata/tdata_n5_mindeg3/X_gt.csv", xGt))
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

    bool readInitguess = true;
    if (readInitguess)
        if (!SomUtils::readCsvInitguess(basepath + "startx.txt", startX))
        {
            ROFL_ERR("Error opening file")
            ROFL_ASSERT(0)
        }

    // RUN RSOM RS
    int srcNodeId = 0;
    SomUtils::VecMatD Rout(n, SomUtils::MatD::Identity(d, d));
    SomUtils::MatD Tout(SomUtils::MatD::Zero(d, n));
    int lastStaircaseStep;
    ROPTLIB::runRsomRS(Prob, startX, srcNodeId, Rout, Tout, lastStaircaseStep); // note: startX is needed (even if random) in ROPTLIB;
    // ROPTLIB namespace is used even if runRsomRS() is not in SampleSomProblem class, nor in "original" ROPTLIB

    std::vector<double> rotErrs(numEdges, 1e+6), translErrs(numEdges, 1e+6);
    ROPTLIB::computeErrorsSingleRsom(edges,
                                     Rout, Tout,
                                     RgtEig, TgtEig,
                                     rotErrs, translErrs);
}