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

    // RUN SSOM RS
    ROPTLIB::Vector randUvec = ProdManiSsom.RandominManifold();
    if (!SomUtils::readCsvInitguess(folderIn + "randUvec.csv", randUvec))
    {
        // matlab/data/ssom_testdata_noisy/harder/tdata_n5_mindeg2_sigma00/ssom_x_start.csv
        ROFL_ERR("Error opening file")
        ROFL_VAR1(folderIn + "randUvec.csv")
        ROFL_ASSERT(0)
    }

    /* Setting up Prob using setters */
    Prob.setRho(5.0);
    Prob.setUsePIM(true);
    Prob.setPimMaxIterations(5000);

    double exectimeGH = 0;
    {
        rofl::ScopedTimer timer("ssom grad");

        startX = -2 * startX;
        double cost = Prob.f(startX);
        ROFL_VAR1(cost)
        startX.Print("startX before grad/hess");

        auto gradOut = Prob.GetDomain()->RandominManifold();

        Prob.Grad(startX, &gradOut);
        // startX.Print("startX after grad");
        // ROFL_VAR1("ssom grad")
        gradOut.Print("gradOut");

        auto hessOut = Prob.GetDomain()->RandominManifold();

        Prob.RieHessianEta(startX, randUvec, &hessOut);
        // startX.Print("startX after grad");
        // ROFL_VAR1("ssom grad")
        hessOut.Print("hessOut");

        exectimeGH = timer.elapsedTimeMs();
    }

    ROFL_VAR1(exectimeGH)

    double exectimePIM = 0;
    {
        rofl::ScopedTimer timer("PIM");

        int staircaseStepIdx = d + 1;
        SomUtils::SomSize somSzNext(staircaseStepIdx, d, n);

        ROPTLIB::SsomProblem ProbNext(somSzNext, Prob.tijs_, Prob.edges_);
        ProbNext.SetDomain(&ProdManiSsom);

        /* Setting up ProbNext using setters */
        ProbNext.setRho(5.0);
        ProbNext.setUsePIM(true);
        ProbNext.setPimMaxIterations(5000);

        ROPTLIB::Stiefel mani1next(staircaseStepIdx, d);
        mani1next.ChooseParamsSet2();
        ROPTLIB::Euclidean mani2next(staircaseStepIdx, n);
        ROPTLIB::ProductManifold ProdManiSsomNext(numoftypes,
                                                  &mani1next, numofmani1, &mani2next, numofmani2, &mani3, numofmani3);
        ROPTLIB::SsomProblem Prob(somSzD, Tijs, edges);

        auto uStart = ProdManiSsomNext.RandominManifold();
        if (!SomUtils::readCsvInitguess(folderIn + "uStartVec.csv", uStart))
        {
            // matlab/data/ssom_testdata_noisy/harder/tdata_n5_mindeg2_sigma00/ssom_x_start.csv
            ROFL_ERR("Error opening file")
            ROFL_VAR1(folderIn + "uStartVec.csv")
            ROFL_ASSERT(0)
        }

        auto uStartSecondIter = ProdManiSsomNext.RandominManifold();
        ProbNext.SetDomain(&ProdManiSsom);
        if (!SomUtils::readCsvInitguess(folderIn + "uStartSecondIterVec.csv", uStartSecondIter))
        {
            // matlab/data/ssom_testdata_noisy/harder/tdata_n5_mindeg2_sigma00/ssom_x_start.csv
            ROFL_ERR("Error opening file")
            ROFL_VAR1(folderIn + "uStartSecondIterVec.csv")
            ROFL_ASSERT(0)
        }

        SomUtils::MatD startXeig(SomUtils::MatD::Zero(d * d * n + d * n + numEdges, 1));
        Prob.RoptToEig(startX, startXeig);
        SomUtils::VecMatD R(SomUtils::VecMatD(n, SomUtils::MatD::Zero(d, d)));
        SomUtils::MatD T(SomUtils::MatD::Zero(d, n));
        SomUtils::MatD Lambdas(SomUtils::MatD::Zero(numEdges, 1));
        Prob.getRotations(startXeig, R);
        Prob.getTranslations(startXeig, T);
        Prob.getScales(startXeig, Lambdas);

        ROPTLIB::Vector Y0out;
        // Y0out = ProbNext.GetDomain()->RandominManifold();
        double lambdaPimOut;

        // double thresh,
        // const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas,
        // Vector &Y0, double &lambdaPimOut,
        // SomUtils::VecMatD &vPimRout, SomUtils::MatD &vPimTout, SomUtils::MatD &vPimLambdasOut,
        // bool armijo = false

        SomUtils::VecMatD vR(n, SomUtils::MatD::Zero(d, d));
        SomUtils::MatD vT(SomUtils::MatD::Zero(d, n));
        SomUtils::MatD vLambdas(SomUtils::MatD::Zero(numEdges, 1));

        ProbNext.ssomPimHessianGenprocEigen(1e-5, R, T, Lambdas, Y0out, lambdaPimOut, vR, vT, vLambdas); //!! catZeroRows() increase is being done inside

        // double lambdaPimMatlab = -384.979;
        // // ProbNext.eigencheckHessianGenproc(lambdaPim, Rnext, vPimR, Tnext, vPimT, LambdasNext, vPimLambdas);
        // // ProbNext.eigencheckHessianGenprocShifted(lambdaPim, Rnext, vPimR, Tnext, vPimT, LambdasNext, vPimLambdas);

        // SomUtils::VecMatD Rnext(n, SomUtils::MatD::Zero(staircaseStepIdx, d));
        // SomUtils::catZeroRow3dArray(R, Rnext);
        // SomUtils::MatD Tnext(SomUtils::MatD::Zero(staircaseStepIdx, n));
        // SomUtils::catZeroRow(T, Tnext);
        // SomUtils::MatD LambdasNext = Lambdas;
        // ProbNext.eigencheckHessianGenproc(lambdaPimMatlab, Rnext, vR, Tnext, vT, LambdasNext, vLambdas);

        SomUtils::VecMatD RnextTgNormStart(n, SomUtils::MatD::Zero(staircaseStepIdx, d));
        SomUtils::MatD TnextTgNormStart(SomUtils::MatD::Zero(staircaseStepIdx, n));
        SomUtils::MatD LambdasNextTgNormStart(SomUtils::MatD::Zero(numEdges, 1));
        SomUtils::MatD uStartEig(SomUtils::MatD::Zero(staircaseStepIdx * d * n + staircaseStepIdx * n + numEdges, 1));
        ProbNext.RoptToEig(uStart, uStartEig);
        ProbNext.getRotations(uStartEig, RnextTgNormStart);
        ProbNext.getTranslations(uStartEig, TnextTgNormStart);
        ProbNext.getScales(uStartEig, LambdasNextTgNormStart);

        SomUtils::VecMatD Rnext2ndTgNormStart(n, SomUtils::MatD::Zero(staircaseStepIdx, d));
        SomUtils::MatD Tnext2ndTgNormStart(SomUtils::MatD::Zero(staircaseStepIdx, n));
        SomUtils::MatD LambdasNext2ndTgNormStart(SomUtils::MatD::Zero(numEdges, 1));
        SomUtils::MatD uStart2ndEig(SomUtils::MatD::Zero(staircaseStepIdx * d * n + staircaseStepIdx * n + numEdges, 1));
        ProbNext.RoptToEig(uStartSecondIter, uStart2ndEig);
        ProbNext.getRotations(uStart2ndEig, Rnext2ndTgNormStart);
        ProbNext.getTranslations(uStart2ndEig, Tnext2ndTgNormStart);
        ProbNext.getScales(uStart2ndEig, LambdasNext2ndTgNormStart);

        // ProbNext.ssomPimHessianGenprocEigenWithStartingPts(1e-5, R, T, Lambdas,
        //                                                    RnextTgNormStart, TnextTgNormStart, LambdasNextTgNormStart,
        //                                                    Rnext2ndTgNormStart, Tnext2ndTgNormStart, LambdasNext2ndTgNormStart,
        //                                                    Y0out, lambdaPimOut, vR, vT, vLambdas, false);

        ProbNext.ssomPimHessianGenprocEigen(1e-5, R, T, Lambdas,
                                            Y0out, lambdaPimOut, vR, vT, vLambdas, false);

        ROFL_VAR1(lambdaPimOut)

        exectimePIM = timer.elapsedTimeMs();
    }

    ROFL_VAR1(exectimePIM)

    return 0;
}
