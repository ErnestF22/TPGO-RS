#include <iostream>
#include <fstream>
#include <set>
#include <filesystem>

#include <rofl/common/param_map.h>
#include <rofl/common/profiler.h>

#include "som_utils.h"

// #include "include/thirdparty/qr_unique_sizeless/include_matlab/lapack.h"
// #include "include/thirdparty/qr_unique_sizeless/include_matlab/lapacke.h"

#include "som_procrustes.h"

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

    if (!SomUtils::readMatlabCsvTijs(folderIn + "/Tijs2.csv", Tijs, d, numEdges))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }
    if (!SomUtils::readMatlabCsvEdges(folderIn + "/edges2.csv", edges))
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

    SomUtils::MatD Xgt = SomUtils::MatD::Random(somSzD.d_ * somSzD.d_ * somSzD.n_ + somSzD.d_ * somSzD.n_, 1);
    SomUtils::readCsvVecEigen(folderIn + "/X2.csv", Xgt);
    std::vector<SomUtils::MatD> Rgt(somSzD.n_, SomUtils::MatD::Identity(somSzD.d_, somSzD.d_));
    SomUtils::MatD RgtVec(SomUtils::MatD::Zero(somSzD.d_ * somSzD.d_ * somSzD.n_, 1));
    RgtVec = Xgt.block(0, 0, somSzD.d_ * somSzD.d_ * somSzD.n_, 1);
    SomUtils::MatD RgtHst = RgtVec.reshaped<Eigen::ColMajor>(somSzD.d_, somSzD.d_ * somSzD.n_);
    SomUtils::unStackH(RgtHst, Rgt, somSzD.d_);

    SomUtils::MatD Tgt = SomUtils::MatD::Random(somSzD.d_ * somSzD.n_, 1);
    Tgt = Xgt.block(somSzD.d_ * somSzD.d_ * somSzD.n_, 0, somSzD.d_ * somSzD.n_, 1);
    Tgt.resize(somSzD.d_, somSzD.n_);
    ROFL_VAR1(Tgt)

    auto Tstart = Tgt; // change this for testing a more realistic case
    sp.setTstart(Tstart);

    SomUtils::MatD A = SomUtils::MatD::Zero(somSzD.d_ * numEdges, somSzD.d_ * somSzD.n_);
    // b=zeros(d*num_edges, 1);
    SomUtils::MatD b = SomUtils::MatD::Zero(somSzD.d_ * numEdges, 1);

    sp.makeAb(Rgt, A, b);

    SomUtils::MatD Amatlab = SomUtils::MatD::Random(somSzD.d_ * numEdges * somSzD.d_ * somSzD.n_, 1);
    SomUtils::MatD bmatlab = SomUtils::MatD::Zero(somSzD.d_ * numEdges, 1);
    SomUtils::readCsvVecEigen(folderIn + "/A_matlab.csv", Amatlab);
    SomUtils::readCsvVecEigen(folderIn + "/b_matlab.csv", bmatlab);
    // Amatlab.resize(somSzD.d_ * numEdges, somSzD.d_ * somSzD.n_);

    std::cout << "A: " << std::endl;
    for (int i = 0; i < Amatlab.rows(); ++i)
        std::cout << A.reshaped<Eigen::ColMajor>(A.rows() * A.cols(), 1)(i, 0) << "\t" << Amatlab(i, 0) << std::endl;
    std::cout << "b: " << std::endl;
    for (int i = 0; i < bmatlab.rows(); ++i)
        std::cout << b.reshaped<Eigen::ColMajor>(b.rows() * b.cols(), 1)(i, 0) << "\t" << bmatlab(i, 0) << std::endl;

    ROFL_VAR2((A.reshaped<Eigen::ColMajor>(A.rows() * A.cols(), 1) - Amatlab).norm(), (b.reshaped<Eigen::ColMajor>(b.rows() * b.cols(), 1) - bmatlab).norm());

    return 0;
}