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

#include <filesystem>

#include "thirdparty/roptlib/problems_SampleSom.h"

namespace fs = std::filesystem;

int main(int argc, char **argv)
{
    std::string filenameCfg;
    std::string folderIn;
    std::set<fs::path> sortedByName;

    int d;
    int numTestsPerInstance;
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
        std::string("../matlab/data/cpp_testdata/"));

    params.getParam<int>("d", d, 3);
    params.getParam<int>("numTestsPerInstance", numTestsPerInstance, 30);
    params.getParam<bool>("readStartingPtFromFile", readStartingPtFromFile, false);
    params.getParam<int>("srcNodeIdx", srcNodeIdx, 0);


    std::cout << "-------\n"
              << std::endl;

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    /*************************End of ROFL params reading**************************/

    for (auto &entry : fs::directory_iterator(folderIn))
        sortedByName.insert(entry.path());

    for (const auto &entry : sortedByName)
    {
        ROFL_VAR1(entry);
        // for (int j = 0; j<numTestsPerInstance; ++j)
        // TODO: repeated tests (e.g., 30) per each test case

        int n;
        SomUtils::readSingleIntCsv(entry.string() + "/n.csv", n);

        int numEdges;

        SomUtils::readSingleIntCsv(entry.string() + "/e.csv", numEdges);

        ROFL_VAR2(n, numEdges); // n = 10, numEdges = 48

        SomUtils::SomSize somSzD(d, d, n);
        SomUtils::MatD Tijs(d, numEdges);
        Eigen::MatrixXi edges(numEdges, 2);

        SomUtils::readMatlabCsvTijs(entry.string() + "/tijs.csv", Tijs, d, numEdges);
        SomUtils::readMatlabCsvEdges(entry.string() + "/edges.csv", edges);

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
        SomUtils::readCsvInitguess(entry.string() + "/X_gt.csv", xGt);
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
            SomUtils::readCsvInitguess("../data/X_initguess_test_single_run_rsom_rs.csv", startX);

        { //rsom RS execution scope
            rofl::ScopedTimer runRsomRStimer ("runRsomRS");

            // RUN RSOM RS
            ROPTLIB::runRsomRS(Prob, startX, srcNodeIdx); // note: startX is needed (even if random) in ROPTLIB;
            // ROPTLIB namespace is used even if runRsomRS() is not in SampleSomProblem class, nor in "original" ROPTLIB

            ROFL_VAR1(runRsomRStimer.elapsedTimeMs())
        }
        
    }

    return 0;
}
