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

void computeErrorsSingleRsom(const Eigen::MatrixXi &edges,
                             const SomUtils::VecMatD &R, const SomUtils::MatD &T,
                             const SomUtils::VecMatD &Rgt, const SomUtils::MatD &Tgt,
                             std::vector<double> &rotErrs, std::vector<double> &translErrs)
{
    // Compute errors
    ROFL_VAR1("Printing R, T out")
    for (auto &m : R)
        ROFL_VAR1(m)
    ROFL_VAR1(T)

    int n = R.size();
    int d = R[0].cols();
    int numEdges = edges.rows();

    for (int e = 0; e < numEdges; ++e)
    {
        int i = edges(e, 0) - 1;
        int j = edges(e, 1) - 1;

        Eigen::Matrix3d ri = R[i].block(0, 0, d, d);
        Eigen::Matrix3d rigt = Rgt[i].block(0, 0, d, d);
        Eigen::Matrix3d rj = R[j].block(0, 0, d, d);
        Eigen::Matrix3d rjgt = Rgt[j].block(0, 0, d, d);
        Eigen::Vector3d ti = T.col(i);
        Eigen::Vector3d tigt = Tgt.col(i);
        Eigen::Vector3d tj = T.col(j);
        Eigen::Vector3d tjgt = Tgt.col(j);

        Eigen::Matrix4d transfI(Eigen::Matrix4d::Identity());
        transfI.block(0, 0, d, d) = ri;
        transfI.block(0, d, d, 1) = ti;
        Eigen::Matrix4d transfJ(Eigen::Matrix4d::Identity());
        transfJ.block(0, 0, d, d) = rj;
        transfJ.block(0, d, d, 1) = tj;

        Eigen::Matrix4d transfIgt(Eigen::Matrix4d::Identity());
        transfIgt.block(0, 0, d, d) = rigt;
        transfIgt.block(0, d, d, 1) = tigt;
        Eigen::Matrix4d transfJgt(Eigen::Matrix4d::Identity());
        transfJgt.block(0, 0, d, d) = rjgt;
        transfJgt.block(0, d, d, 1) = tjgt;

        Eigen::MatrixXd p(Eigen::MatrixXd::Identity(d + 1, d + 1));
        SomUtils::computeRelativePose(transfI, transfJ, p);

        Eigen::MatrixXd pGt(Eigen::MatrixXd::Identity(d + 1, d + 1));
        SomUtils::computeRelativePose(transfIgt, transfJgt, pGt);

        Eigen::Matrix3d pR = p.block(0, 0, d, d);
        Eigen::Matrix3d pRgt = pGt.block(0, 0, d, d);

        double rotDistEdge = SomUtils::rotDistSingle(pR, pRgt);
        ROFL_VAR2(e, rotDistEdge);

        Eigen::Vector3d pT = p.block(0, d, d, 1);
        Eigen::Vector3d pTgt = pGt.block(0, d, d, 1);

        double translDistEdge = SomUtils::translErr(ti, tigt);
        ROFL_VAR2(ti.transpose(), tigt.transpose());
        ROFL_VAR2(e, translDistEdge);

        rotErrs[e] = rotDistEdge;
        translErrs[e] = translDistEdge;
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
    for (int i = 0; i< Amatlab.rows(); ++i)
        std::cout << A.reshaped<Eigen::ColMajor>(A.rows() * A.cols(), 1)(i,0) << "\t" << Amatlab(i,0) << std::endl;
    std::cout << "b: " << std::endl;
    for (int i = 0; i< bmatlab.rows(); ++i)
        std::cout << b.reshaped<Eigen::ColMajor>(b.rows() * b.cols(), 1)(i,0) << "\t" << bmatlab(i,0) << std::endl;

    ROFL_VAR2((A.reshaped<Eigen::ColMajor>(A.rows() * A.cols(), 1) - Amatlab).norm(), (b.reshaped<Eigen::ColMajor>(b.rows() * b.cols(), 1) - bmatlab).norm());

    return 0;
}