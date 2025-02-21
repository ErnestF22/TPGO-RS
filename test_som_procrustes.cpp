#include <iostream>
#include <fstream>
#include <set>
#include <filesystem>

#include <rofl/common/param_map.h>
#include <rofl/common/profiler.h>

// #include "include/thirdparty/qr_unique_sizeless/include_matlab/lapack.h"
// #include "include/thirdparty/qr_unique_sizeless/include_matlab/lapacke.h"

#include "som_procrustes.h"

namespace fs = std::filesystem;

void readCsvVecEigen(const std::string &filenameIn, Eigen::MatrixXd &out)
{
    std::string line;
    std::ifstream fileIn(filenameIn);
    if (!fileIn.is_open())
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }
    int ctr = 0;
    while (std::getline(fileIn, line))
    {
        // ROFL_VAR1(line);
        double val = std::stod(line);
        out(ctr, 0) = val;
        ctr++;
    }
}

int main(int argc, char **argv)
{

    std::string filenameCfg;
    std::string folderIn;
    std::set<fs::path> sortedByName;

    int d;
    int numTestsPerInstance;
    bool readStartingPtFromFile;
    int srcNodeIdx;
    std::string resultsBasePath;
    std::string fileTstart;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>(
        "in", folderIn,
        std::string("../matlab/data/cpp_testdata/tdata_n5_mindeg3"));

    params.getParam<int>("d", d, 3);
    params.getParam<int>("numTestsPerInstance", numTestsPerInstance, 30);
    params.getParam<bool>("readStartingPtFromFile", readStartingPtFromFile, false);
    params.getParam<std::string>("resultsBasePath", resultsBasePath, "../results_procrustes/");
    params.getParam<int>("srcNodeIdx", srcNodeIdx, 0);
    params.getParam<std::string>("fileTstart", fileTstart, "../data/procrustes_t_start.csv");


    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    std::cout << "-------\n"
              << std::endl;

    /*************************End of ROFL params reading**************************/

    // ROFL_VAR1(entry);
    // for (int j = 0; j<numTestsPerInstance; ++j)
    // TODO: repeated tests (e.g., 30) per each test case

    int n;
    if (!SomUtils::readSingleIntCsv(folderIn + "/n.csv", n))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }

    int numEdges;

    if (!SomUtils::readSingleIntCsv(folderIn + "/e.csv", numEdges))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }

    ROFL_VAR2(n, numEdges); // n = 5, numEdges = 18

    SomUtils::SomSize somSzD(d, d, n);
    SomUtils::MatD Tijs(d, numEdges);
    Eigen::MatrixXi edges(numEdges, 2);

    if (!SomUtils::readMatlabCsvTijs(folderIn + "/tijs.csv", Tijs, d, numEdges))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }
    if (!SomUtils::readMatlabCsvEdges(folderIn + "/edges.csv", edges))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }

    ROFL_VAR1(Tijs)
    ROFL_VAR1(edges)

    // gt
    //  Read GT from csv
    //  ROPTLIB::Vector xGt = ProdMani.RandominManifold();
    //  if (!SomUtils::readCsvInitguess(entry.string() + "/Xgt.csv", xGt))
    //  {
    //      ROFL_ERR("Error opening file")
    //      ROFL_ASSERT(0)
    //  }
    //  xGt.Print("xGt");
    //  // ROPT to Eig (GT)
    //  SomUtils::MatD XgtVecEig(SomUtils::MatD::Zero(d * d * n + d * n, 1));
    //  SomUtils::VecMatD RgtEig(n, SomUtils::MatD::Zero(d, d));
    //  SomUtils::MatD TgtEig(SomUtils::MatD::Zero(d, n));
    //  Prob.RoptToEig(xGt, XgtVecEig);
    //  Prob.getRotations(XgtVecEig, RgtEig);
    //  Prob.getTranslations(XgtVecEig, TgtEig);
    //  Prob.setGt(RgtEig, TgtEig);

    // method initialization and run
    SomProcrustes sp(somSzD, Tijs, edges);

    SomUtils::MatD Xstart = SomUtils::MatD::Random(somSzD.d_ * somSzD.d_ * somSzD.n_ + somSzD.d_ * somSzD.n_, 1);
    readCsvVecEigen(folderIn + "/X_gt.csv", Xstart);
    SomUtils::MatD Tstart = SomUtils::MatD::Random(somSzD.d_ * somSzD.n_, 1);
    Tstart = Xstart.block(somSzD.d_ * somSzD.d_ * somSzD.n_, 0, somSzD.d_ * somSzD.n_, 1);
    Tstart.resize(somSzD.d_, somSzD.n_);
    ROFL_VAR1(Tstart)
    sp.setTcurr(Tstart);

    sp.run();

    SomUtils::MatD Tout = sp.getTout();
    SomUtils::VecMatD Rout = sp.getRout();

    for (int i = 0; i < n; ++i)
    {
        ROFL_VAR2(i, Rout[i]);
    }
    ROFL_VAR1(Tout);

    return 0;
}