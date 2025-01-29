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
    std::string resultsBasePath;

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
    params.getParam<std::string>("resultsBasePath", resultsBasePath, "../results/");
    params.getParam<int>("srcNodeIdx", srcNodeIdx, 0);

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    std::cout << "-------\n"
              << std::endl;

    /*************************End of ROFL params reading**************************/

    for (auto &entry : fs::directory_iterator(folderIn))
        sortedByName.insert(entry.path());

    std::vector<std::vector<std::vector<double>>> rotErrsAll, translErrsAll;
    std::vector<std::vector<double>> execTimesAll;
    std::vector<std::vector<double>> staircaseStepOutIdxAll;

    for (const auto &entry : sortedByName)
    {
        std::vector<std::vector<double>> rotErrs(numTestsPerInstance), translErrs(numTestsPerInstance);
        std::vector<double> execTimes(numTestsPerInstance), staircaseStepOutIdx(numTestsPerInstance);

        // ROFL_VAR1(entry);
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

        if (!fs::exists(resultsBasePath))
            fs::create_directory(resultsBasePath);

        int pos = entry.string().find("mindeg");
        std::string mindegStr = entry.string().substr(pos + 6, 1); // mindeg has 6 characters
        ROFL_VAR1(mindegStr);
        int mindeg = boost::lexical_cast<int, std::string>(mindegStr);
        ROFL_VAR1(mindeg)

        std::string folderAppendName =
            "n" + boost::lexical_cast<std::string, int>(n) +
            "_mindeg" + boost::lexical_cast<std::string, int>(mindeg) +
            "_noise00"; // TODO: make noise param reading automated
        std::string folderAppendNameStamped = SomUtils::generateStampedString(folderAppendName, "/");
        fs::create_directory(resultsBasePath + folderAppendNameStamped);

        // 1
        std::string rotErrsFilename = resultsBasePath + folderAppendNameStamped + folderAppendName + "_rot_errors.txt";
        std::ofstream rotErrsOfstream;
        rotErrsOfstream.open(rotErrsFilename);
        std::string translErrsFilename = resultsBasePath + folderAppendNameStamped + folderAppendName + "_transl_errors.txt";
        std::ofstream translErrsOfstream;
        translErrsOfstream.open(translErrsFilename);
        std::string execTimesFilename = resultsBasePath + folderAppendNameStamped + folderAppendName + "_exec_times.txt";
        std::ofstream execTimesOfstream;
        execTimesOfstream.open(execTimesFilename);
        std::string staircaseStepOutIdxFilename = resultsBasePath + folderAppendNameStamped + folderAppendName + "_last_rs_step.txt";
        std::ofstream staircaseStepOutIdxOfstream;
        staircaseStepOutIdxOfstream.open(staircaseStepOutIdxFilename);
        // means
        std::string rotErrsMeanFilename = resultsBasePath + folderAppendNameStamped + folderAppendName + "_rot_errors_mean.txt";
        std::ofstream rotErrsMeanOfstream;
        rotErrsMeanOfstream.open(rotErrsMeanFilename);
        std::string translErrsMeanFilename = resultsBasePath + folderAppendNameStamped + folderAppendName + "_transl_errors_mean.txt";
        std::ofstream translErrsMeanOfstream;
        translErrsMeanOfstream.open(translErrsMeanFilename);
        std::string execTimesMeanFilename = resultsBasePath + folderAppendNameStamped + folderAppendName + "_exec_times_mean.txt";
        std::ofstream execTimesMeanOfstream;
        execTimesMeanOfstream.open(execTimesMeanFilename);

        for (int testjd = 0; testjd < numTestsPerInstance; ++testjd)
        {
            ROFL_VAR3(entry, testjd, "start");

            // Generate startX (random)
            ROPTLIB::Vector startX = ProdMani.RandominManifold();

            if (readStartingPtFromFile)
                SomUtils::readCsvInitguess("../data/X_initguess_test_single_run_rsom_rs.csv", startX);

            startX.Print("startX");

            { // rsom RS execution scope
                rofl::ScopedTimer runRsomRStimer("runRsomRS");

                // RUN RSOM RS
                SomUtils::VecMatD Rout(n, SomUtils::MatD::Identity(d, d));
                SomUtils::MatD Tout(SomUtils::MatD::Zero(d, n));
                int staircaseStepIdxOutIJ;
                ROPTLIB::runRsomRS(Prob, startX, srcNodeIdx, Rout, Tout, staircaseStepIdxOutIJ); // note: startX is needed (even if random) in ROPTLIB;
                // ROPTLIB namespace is used even if runRsomRS() is not in SampleSomProblem class, nor in "original" ROPTLIB

                std::vector<double> rotErrsTestjd(numEdges, 1e+6), translErrsTestjd(numEdges, 1e+6);
                ROPTLIB::computeErrorsSingleRsom(edges,
                                                 Rout, Tout,
                                                 RgtEig, TgtEig,
                                                 rotErrsTestjd, translErrsTestjd);

                double execTimeIJ = runRsomRStimer.elapsedTimeMs();
                ROFL_VAR3(entry, testjd, execTimeIJ)

                rotErrs[testjd].resize(numEdges);
                translErrs[testjd].resize(numEdges);

                rotErrs[testjd] = rotErrsTestjd;
                translErrs[testjd] = translErrsTestjd;

                execTimes[testjd] = execTimeIJ;
                staircaseStepOutIdx[testjd] = staircaseStepIdxOutIJ;
            } // end of rsom RS execution scope
            // break; //testjd loop
        }
        // break; //entries/instances loop
        for (int i = 0; i < numTestsPerInstance; ++i)
        {
            for (int j = 0; j < numEdges; ++j)
            {
                // output 2
                // ROFL_VAR4(entry, j, rotErrs[i][j], translErrs[i][j]);
                rotErrsOfstream << "i " + std::to_string(i) + " j " << std::to_string(j) << std::endl;
                rotErrsOfstream << rotErrs[i][j] << std::endl;
                translErrsOfstream << "i " + std::to_string(i) + " j " << std::to_string(j) << std::endl;
                translErrsOfstream << translErrs[i][j] << std::endl;
            }
            execTimesOfstream << "i " + std::to_string(i) << std::endl;
            execTimesOfstream << execTimes[i] << std::endl;
            staircaseStepOutIdxOfstream << "i " + std::to_string(i) << std::endl;
            staircaseStepOutIdxOfstream << staircaseStepOutIdx[i] << std::endl;


            // Finding mean error of current instance-testjd pair
            double rotMeanErr = std::accumulate(rotErrs[i].begin(), rotErrs[i].end(), 0) / rotErrs[i].size();
            double translMeanErr = std::accumulate(translErrs[i].begin(), translErrs[i].end(), 0) / translErrs[i].size();
            rotErrsMeanOfstream << "i " + std::to_string(i) << std::endl;
            rotErrsMeanOfstream << rotMeanErr << std::endl;
            translErrsMeanOfstream << "i " + std::to_string(i) << std::endl;
            translErrsMeanOfstream << translMeanErr << std::endl;
        }

        double exectimeMean = std::accumulate(execTimes.begin(), execTimes.end(), 0) / execTimes.size();
        // execTimesMeanOfstream << "i " + std::to_string(i) << std::endl;
        execTimesMeanOfstream << exectimeMean << std::endl;

        rotErrsAll.push_back(rotErrs);
        translErrsAll.push_back(translErrs);
        execTimesAll.push_back(execTimes);

        // output 3
        rotErrsOfstream.close();
        translErrsOfstream.close();
        execTimesOfstream.close();
        rotErrsMeanOfstream.close();
        translErrsMeanOfstream.close();
        execTimesMeanOfstream.close();
    }

    // for (int i = 0; i < sortedByName.size(); ++i)
    // {
    //     for (int j = 0; j < numTestsPerInstance; ++j)
    //         for (int k = 0; k < rotErrsAll[i][j].size(); ++k)
    //             ROFL_VAR5(i, j, k, rotErrsAll[i][j][k], translErrsAll[i][j][k]);
    // }

    return 0;
}
