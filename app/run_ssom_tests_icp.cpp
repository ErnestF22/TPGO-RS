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

#include "thirdparty/roptlib/problems_Ssom.h"

#include "som_utils.h"

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
    double rho;

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
    params.getParam<int>("numTestsPerInstance", numTestsPerInstance, 2);
    params.getParam<bool>("readStartingPtFromFile", readStartingPtFromFile, true);
    params.getParam<std::string>("resultsBasePath", resultsBasePath, "../results_ssom_icp/");
    params.getParam<int>("srcNodeIdx", srcNodeIdx, 0);
    params.getParam<double>("rho", rho, 0.0);

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    std::cout << "-------\n"
              << std::endl;

    /*************************End of ROFL params reading**************************/

    for (auto &entry : fs::directory_iterator(folderIn))
        sortedByName.insert(entry.path());

    std::vector<std::vector<std::vector<double>>> rotErrsAll, translErrsAll, LambdasMeanErrsAll, LambdasMaxErrsAll;
    std::vector<std::vector<double>> execTimesAll;

    // declaring ofstreams
    std::ofstream rotErrsOfstream;    
    std::ofstream translErrsOfstream;
    std::ofstream rotErrsMeanOfstream;
    std::ofstream translErrsMeanOfstream;
    std::ofstream execTimesOfstream;
    std::ofstream execTimesMeanOfstream;
    std::ofstream LambdasMeanErrsOfstream;
    std::ofstream LambdasMaxErrsOfstream;
    std::ofstream LambdasMeanErrsMeanOfstream;
    std::ofstream LambdasMaxErrsMeanOfstream;

    std::string folderAppendNameStamped = SomUtils::generateStampedString("", "");

    int numInstances = sortedByName.size();

    int inst = 0;
    for (const auto &entry : sortedByName)
    {
        // ROFL_VAR1(entry);
        // for (int j = 0; j<numTestsPerInstance; ++j)
        // TODO: repeated tests (e.g., 30) per each test case

        int n;
        if (!SomUtils::readSingleIntCsv(entry.string() + "/n.csv", n))
        {
            ROFL_ERR("Error opening file")
            ROFL_ASSERT(0)
        }

        int numEdges;

        if (!SomUtils::readSingleIntCsv(entry.string() + "/e.csv", numEdges))
        {
            ROFL_ERR("Error opening file")
            ROFL_ASSERT(0)
        }

        ROFL_VAR2(n, numEdges); // n = 10, numEdges = 48

        SomUtils::SomSize somSzD(d, d, n);
        SomUtils::MatD Tijs(d, numEdges);
        Eigen::MatrixXi edges(numEdges, 2);

        if (!SomUtils::readMatlabCsvTijs(entry.string() + "/tijs.csv", Tijs, d, numEdges))
        {
            ROFL_ERR("Error opening file")
            ROFL_ASSERT(0)
        }
        if (!SomUtils::readMatlabCsvEdges(entry.string() + "/edges.csv", edges))
        {
            ROFL_ERR("Error opening file")
            ROFL_ASSERT(0)
        }

        ROFL_VAR1(Tijs)
        ROFL_VAR1(edges)

        int e = numEdges;

        // problem d x d x n
        integer numoftypes = 3; // 2 i.e. (3D) Stiefel + Euclidean
        integer numofmani1 = n; // num of Stiefel manifolds
        integer numofmani2 = 1;
        integer numofmani3 = 1;
        ROPTLIB::Stiefel mani1(d, d);
        mani1.ChooseParamsSet2();
        ROPTLIB::Euclidean mani2(d, n);
        ROPTLIB::Euclidean mani3(e);
        ROPTLIB::ProductManifold ProdMani(numoftypes, &mani1, numofmani1, &mani2, numofmani2, &mani3, numofmani3);
        ROPTLIB::SsomProblem Prob(somSzD, Tijs, edges);

        // Read GT from csv
        ROPTLIB::Vector xGt = ProdMani.RandominManifold();
        if (!SomUtils::readCsvInitguess(entry.string() + "/Xgt.csv", xGt))
        {
            ROFL_ERR("Error opening file")
            ROFL_ASSERT(0)
        }
        xGt.Print("xGt");
        // ROPT to Eig (GT)
        SomUtils::MatD XgtVecEig(SomUtils::MatD::Zero(d * d * n + d * n + e, 1));
        SomUtils::VecMatD RgtEig(n, SomUtils::MatD::Zero(d, d));
        SomUtils::MatD TgtEig(SomUtils::MatD::Zero(d, n));
        SomUtils::MatD LambdasGtEig(SomUtils::MatD::Zero(e, 1));
        Prob.RoptToEig(xGt, XgtVecEig);
        Prob.getRotations(XgtVecEig, RgtEig);
        Prob.getTranslations(XgtVecEig, TgtEig);
        Prob.getScales(XgtVecEig, LambdasGtEig);
        Prob.setGt(RgtEig, TgtEig, LambdasGtEig);
        Prob.setRho(rho);

        // Set the domain of the problem to be the product of Stiefel manifolds
        Prob.SetDomain(&ProdMani);

        // Set Prob params
        Prob.SetUseGrad(true);
        Prob.SetUseHess(true);

        ROFL_VAR1("Printing R, T gt")
        for (auto &m : RgtEig)
            ROFL_VAR1(m)
        ROFL_VAR1(TgtEig)
        ROFL_VAR1(LambdasGtEig)
        // // problem nrs x d x n
        // ROPTLIB::Stiefel mani1nrs(nrs, d);
        // mani1.ChooseParamsSet2();
        // ROPTLIB::Euclidean mani2nrs(nrs, n);
        // ROPTLIB::ProductManifold ProdManiNrs(numoftypes, &mani1nrs, numofmani1, &mani2nrs, numofmani2);
        // // Read from csv
        // ROPTLIB::Vector xManoptOut = ProdManiNrs.RandominManifold();
        // SomUtils::readCsvInitguess("../data/recov_x_manopt_out.csv", xManoptOut);
        // xManoptOut.Print("xManoptOut");

        if (!fs::exists(resultsBasePath))
            fs::create_directory(resultsBasePath);

        int pos = entry.string().find("mindeg");
        std::string mindegStr = entry.string().substr(pos + 6, 1); // mindeg has 6 characters
        ROFL_VAR1(mindegStr);
        int mindeg = boost::lexical_cast<int, std::string>(mindegStr);
        ROFL_VAR1(mindeg)

        int pos2 = entry.string().find("sigma");
        std::string sigmaStr = entry.string().substr(pos2 + 5, 3); // sigma has 5 characters
        ROFL_VAR1(sigmaStr);

        ROFL_VAR2(n, mindeg)

        std::string folderAppendName =
            "n" + boost::lexical_cast<std::string, int>(n) +
            "_mindeg" + boost::lexical_cast<std::string, int>(mindeg) +
            "_sigma" + sigmaStr; // TODO: make noise param reading automated
        fs::create_directory(resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/");

        // 1
        std::string rotErrsFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_rot_errors.txt";

        if (rotErrsOfstream.is_open())
        {
            rotErrsOfstream.close();
            rotErrsOfstream.clear(); // clear flags
        }
        rotErrsOfstream.open(rotErrsFilename);
        if (!rotErrsOfstream)
        {
            ROFL_ERR("Error opening output file")
            ROFL_ASSERT(0)
        }
        std::string translErrsFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_transl_errors.txt";
        if (translErrsOfstream.is_open())
        {
            translErrsOfstream.close();
            translErrsOfstream.clear(); // clear flags
        }
        translErrsOfstream.open(translErrsFilename);
        if (!translErrsOfstream)
        {
            ROFL_ERR("Error opening output file")
            ROFL_ASSERT(0)
        }
        std::string LambdasMeanErrsFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_lambda_mean_errors.txt";
        if (LambdasMeanErrsOfstream.is_open())
        {
            LambdasMeanErrsOfstream.close();
            LambdasMeanErrsOfstream.clear(); // clear flags
        }
        LambdasMeanErrsOfstream.open(LambdasMeanErrsFilename);
        if (!LambdasMeanErrsOfstream)
        {
            ROFL_ERR("Error opening output file")
            ROFL_ASSERT(0)
        }
        std::string LambdasMaxErrsFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_lambda_max_errors.txt";
        if (LambdasMaxErrsOfstream.is_open())
        {
            LambdasMaxErrsOfstream.close();
            LambdasMaxErrsOfstream.clear(); // clear flags
        }
        LambdasMaxErrsOfstream.open(LambdasMaxErrsFilename);
        if (!LambdasMaxErrsOfstream)
        {
            ROFL_ERR("Error opening output file")
            ROFL_ASSERT(0)
        }
        std::string execTimesFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_exec_times.txt";
        if (execTimesOfstream.is_open())
        {
            execTimesOfstream.close();
            execTimesOfstream.clear(); // clear flags
        }
        execTimesOfstream.open(execTimesFilename);
        if (!execTimesOfstream)
        {
            ROFL_ERR("Error opening output file")
            ROFL_ASSERT(0)
        }

        // means
        std::string rotErrsMeanFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_rot_errors_mean.txt";
        if (rotErrsMeanOfstream.is_open())
        {
            rotErrsMeanOfstream.close();
            rotErrsMeanOfstream.clear(); // clear flags
        }
        rotErrsMeanOfstream.open(rotErrsMeanFilename);
        if (!rotErrsMeanOfstream)
        {
            ROFL_ERR("Error opening output file")
            ROFL_ASSERT(0)
        }
        std::string translErrsMeanFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_transl_errors_mean.txt";
        if (translErrsMeanOfstream.is_open())
        {
            translErrsMeanOfstream.close();
            translErrsMeanOfstream.clear(); // clear flags
        }
        translErrsMeanOfstream.open(translErrsMeanFilename);
        if (!translErrsMeanOfstream)
        {
            ROFL_ERR("Error opening output file")
            ROFL_ASSERT(0)
        }
        std::string LambdasMeanErrsMeanFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_lambda_mean_errors_mean.txt";
        if (LambdasMeanErrsMeanOfstream.is_open())
        {
            LambdasMeanErrsMeanOfstream.close();
            LambdasMeanErrsMeanOfstream.clear(); // clear flags
        }
        LambdasMeanErrsMeanOfstream.open(LambdasMeanErrsMeanFilename);
        if (!LambdasMeanErrsMeanOfstream)
        {
            ROFL_ERR("Error opening output file")
            ROFL_ASSERT(0)
        }
        std::string LambdasMaxErrsMeanFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_lambda_max_errors_mean.txt";
        if (LambdasMaxErrsMeanOfstream.is_open())
        {
            LambdasMaxErrsMeanOfstream.close();
            LambdasMaxErrsMeanOfstream.clear(); // clear flags
        }
        LambdasMaxErrsMeanOfstream.open(LambdasMaxErrsMeanFilename);
        if (!LambdasMaxErrsMeanOfstream)
        {
            ROFL_ERR("Error opening output file")
            ROFL_ASSERT(0)
        }
        std::string execTimesMeanFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_exec_times_mean.txt";
        if (execTimesMeanOfstream.is_open())
        {
            execTimesMeanOfstream.close();
            execTimesMeanOfstream.clear(); // clear flags
        }
        execTimesMeanOfstream.open(execTimesMeanFilename);
        if (!execTimesMeanOfstream)
        {
            ROFL_ERR("Error opening output file")
            ROFL_ASSERT(0)
        }

        // Declare and init error metrics (for each instances)
        double rotMeanErr = 1e+6, translMeanErr = 1e+6, execTimeMean = 1e+6, LambdasMeanErr = 1e+6, LambdasMaxErr = 1e+6;
        std::vector<std::vector<double>> rotErrs(numTestsPerInstance), translErrs(numTestsPerInstance), LambdasMeanErrs(numTestsPerInstance), LambdasMaxErrs(numTestsPerInstance);
        std::vector<double> execTimes(numTestsPerInstance);

        for (int testjd = 0; testjd < numTestsPerInstance; ++testjd)
        {
            ROFL_VAR3(entry, testjd, "start");

            // Generate startX (random)
            ROPTLIB::Vector startX = ProdMani.RandominManifold();
            // startX.Initialization(numoftypes, &mani1, numofmani1, &mani2, numofmani2);

            if (readStartingPtFromFile)
                if (!SomUtils::readCsvInitguess(entry.string() + "/ssom_x_start.csv", startX))
                {
                    ROFL_ERR("Error opening file")
                    ROFL_ASSERT(0)
                }

            startX.Print("startX");

            { // Ssom RS execution scope
                rofl::ScopedTimer runSsomICPtimer("runSsomICP");

                // RUN SSOM RS
                SomUtils::VecMatD Rout(n, SomUtils::MatD::Identity(d, d));
                SomUtils::MatD Tout(SomUtils::MatD::Zero(d, n));
                SomUtils::MatD LambdasOut(SomUtils::MatD::Zero(e, 1));
                double costOut = ROPTLIB::runSsomICP(Prob, startX, srcNodeIdx, Rout, Tout, LambdasOut); // note: startX is needed (even if random) in ROPTLIB;
                // ROPTLIB namespace is used even if runSsomRS() is not in SampleSomProblem class, nor in "original" ROPTLIB

                // costOut = 1.0f;
                if (!SomUtils::isEqualDoubles(costOut, 0.0f) && sigmaStr == "00" && testjd == 0)
                {
                    auto fstr = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_j" + std::to_string(testjd);
                    ROFL_VAR1(costOut)
                    std::ofstream tijsofs(fstr + "_tijs.txt");
                    tijsofs << Tijs;
                    std::ofstream edgesofs(fstr + "_edges.txt");
                    edgesofs << edges;
                    std::ofstream startxofs(fstr + "_startx.txt");
                    SomUtils::MatD startXeig(SomUtils::MatD::Zero(d * d * n + d * n + e, 1));
                    Prob.RoptToEig(startX, startXeig);
                    startxofs << startXeig;
                }

                std::vector<double> rotErrsTestjd(numEdges, 1e+6), translErrsTestjd(numEdges, 1e+6), LambdasMeanErrsTestjd(numEdges, 1e+6), LambdasMaxErrsTestjd(numEdges, 1e+6);
                SomUtils::computeErrorsSingleSsom(edges,
                                                  Rout, Tout, LambdasOut,
                                                  RgtEig, TgtEig, LambdasGtEig,
                                                  rotErrsTestjd, translErrsTestjd, LambdasMeanErrsTestjd, LambdasMaxErrsTestjd);

                double execTimeIJ = runSsomICPtimer.elapsedTimeMs();
                ROFL_VAR3(entry, testjd, execTimeIJ)

                rotErrs[testjd].resize(numEdges, 1e+6);
                translErrs[testjd].resize(numEdges, 1e+6);
                LambdasMeanErrs[testjd].resize(numEdges, 1e+6);
                LambdasMaxErrs[testjd].resize(numEdges, 1e+6);

                rotErrs[testjd] = rotErrsTestjd;
                translErrs[testjd] = translErrsTestjd;
                LambdasMeanErrs[testjd] = LambdasMeanErrsTestjd;
                LambdasMaxErrs[testjd] = LambdasMaxErrsTestjd;

                execTimes[testjd] = execTimeIJ;

                for (int k = 0; k < numEdges; ++k)
                {
                    // output 2
                    // ROFL_VAR4(entry, j, rotErrs[i][j], translErrs[i][j]);
                    rotErrsOfstream << "testjd " + std::to_string(testjd) + " k " << std::to_string(k) << std::endl;
                    rotErrsOfstream << rotErrs[testjd][k] << std::endl;
                    translErrsOfstream << "testjd " + std::to_string(testjd) + " k " << std::to_string(k) << std::endl;
                    translErrsOfstream << translErrs[testjd][k] << std::endl;
                    LambdasMeanErrsOfstream << "testjd " + std::to_string(testjd) + " k " << std::to_string(k) << std::endl;
                    LambdasMeanErrsOfstream << LambdasMeanErrs[testjd][k] << std::endl;
                    LambdasMaxErrsOfstream << "testjd " + std::to_string(testjd) + " k " << std::to_string(k) << std::endl;
                    LambdasMaxErrsOfstream << LambdasMaxErrs[testjd][k] << std::endl;
                    ROFL_VAR4(rotErrs[testjd][k], translErrs[testjd][k], LambdasMeanErrs[testjd][k], LambdasMaxErrs[testjd][k])
                }
                execTimesOfstream << "j " + std::to_string(testjd) << std::endl;
                execTimesOfstream << execTimes[testjd] << std::endl;

                // Finding mean error of current instance-testjd pair
                rotMeanErr = SomUtils::stlVecDoublesMean(rotErrs[testjd]);
                translMeanErr = SomUtils::stlVecDoublesMean(translErrs[testjd]);
                LambdasMeanErr = SomUtils::stlVecDoublesMean(LambdasMeanErrs[testjd]);
                LambdasMaxErr = SomUtils::stlVecDoublesMean(LambdasMaxErrs[testjd]);
                rotErrsMeanOfstream << "j " + std::to_string(testjd) << std::endl;
                rotErrsMeanOfstream << rotMeanErr << std::endl;
                translErrsMeanOfstream << "j " + std::to_string(testjd) << std::endl;
                translErrsMeanOfstream << translMeanErr << std::endl;
                LambdasMeanErrsOfstream << "j " + std::to_string(testjd) << std::endl;
                LambdasMeanErrsOfstream << LambdasMeanErr << std::endl;
                LambdasMaxErrsOfstream << "j " + std::to_string(testjd) << std::endl;
                LambdasMaxErrsOfstream << LambdasMaxErr << std::endl;
                ROFL_VAR4(rotMeanErr, translMeanErr, LambdasMeanErr, LambdasMaxErr)

                execTimeMean = SomUtils::stlVecDoublesMean(execTimes);
            } // end of rsom RS execution scope
            // break; //testjd loop
            // startX.Delete();
        } // end of for testjd = 0 : numTestsPerInstance
        // break; //entries/instances loop
        for (int j = 0; j < numTestsPerInstance; ++j)
        {
            for (int k = 0; k < numEdges; ++k)
            {
                // output 2
                // ROFL_VAR4(entry, j, rotErrs[i][j], translErrs[i][j]);
                ROFL_VAR4(rotErrs[j][k], translErrs[j][k], LambdasMeanErrs[j][k], LambdasMaxErrs[j][k])
            }
            ROFL_VAR4(rotMeanErr, translMeanErr, LambdasMeanErr, LambdasMaxErr)
        }

        // execTimesMeanOfstream << "i " + std::to_string(i) << std::endl;
        execTimesMeanOfstream << execTimeMean << std::endl;

        rotErrsAll.push_back(rotErrs);
        translErrsAll.push_back(translErrs);
        LambdasMeanErrsAll.push_back(LambdasMeanErrs);
        LambdasMaxErrsAll.push_back(LambdasMaxErrs);
        execTimesAll.push_back(execTimes);

        inst++; // current instance idx
    }

    rotErrsOfstream.close();
    translErrsOfstream.close();
    LambdasMeanErrsOfstream.close();
    LambdasMaxErrsOfstream.close();
    execTimesOfstream.close();
    rotErrsMeanOfstream.close();
    translErrsMeanOfstream.close();
    LambdasMeanErrsMeanOfstream.close();
    LambdasMaxErrsMeanOfstream.close();
    execTimesMeanOfstream.close();

    for (int i = 0; i < numInstances; ++i) // i already declared
    {
        for (int j = 0; j < numTestsPerInstance; ++j)
        {
            if (!rotErrsAll.empty())
            {
                for (int k = 0; k < rotErrsAll[i][j].size(); ++k)
                    ROFL_VAR7(i, j, k, rotErrsAll[i][j][k], translErrsAll[i][j][k], LambdasMeanErrsAll[i][j][k], LambdasMaxErrsAll[i][j][k]);
            }
            else
                continue;

            ROFL_VAR5(SomUtils::stlVecDoublesMean(rotErrsAll[i][j]), SomUtils::stlVecDoublesMean(translErrsAll[i][j]), SomUtils::stlVecDoublesMean(LambdasMeanErrsAll[i][j]), SomUtils::stlVecDoublesMean(LambdasMaxErrsAll[i][j]),    
                      SomUtils::stlVecDoublesMean(execTimesAll[i]));
        }
    }

    return 0;
}
