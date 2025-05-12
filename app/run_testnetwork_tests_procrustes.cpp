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

#include "som_procrustes.h"
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

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>(
        "in", folderIn,
        std::string("../matlab/data/cpp_testdata_noisy/"));

    params.getParam<int>("d", d, 3);
    params.getParam<int>("numTestsPerInstance", numTestsPerInstance, 30);
    params.getParam<bool>("readStartingPtFromFile", readStartingPtFromFile, true);
    params.getParam<std::string>("resultsBasePath", resultsBasePath, "../results_procrustes/");
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

    // declaring ofstreams
    std::ofstream rotErrsOfstream;
    std::ofstream translErrsOfstream;
    std::ofstream execTimesOfstream;
    std::ofstream rotErrsMeanOfstream;
    std::ofstream translErrsMeanOfstream;
    std::ofstream execTimesMeanOfstream;

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

        if (!fs::exists(resultsBasePath))
            fs::create_directory(resultsBasePath);

        int pos = entry.string().find("mindeg");
        std::string mindegStr = entry.string().substr(pos + 6, 1); // mindeg has 6 characters
        ROFL_VAR1(mindegStr);
        int mindeg = boost::lexical_cast<int, std::string>(mindegStr);
        ROFL_VAR1(mindeg)

        int pos2 = entry.string().find("sigma");
        std::string sigmaStr = entry.string().substr(pos2 + 5, 2); // sigma has 5 characters
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
        double rotMeanErr = 1e+6, translMeanErr = 1e+6, execTimeMean = 1e+6;
        std::vector<std::vector<double>> rotErrs(numTestsPerInstance), translErrs(numTestsPerInstance);
        std::vector<double> execTimes(numTestsPerInstance);

        for (int testjd = 0; testjd < numTestsPerInstance; ++testjd)
        {
            ROFL_VAR3(entry, testjd, "start");

            // Generate startX (random)
            SomUtils::MatD startX = SomUtils::MatD::Random(somSzD.d_ * somSzD.d_ * somSzD.n_ + somSzD.d_ * somSzD.n_, 1);

            if (readStartingPtFromFile)
                if (!SomUtils::readCsvVecEigen(entry.string() + "/startX.csv", startX))
                {
                    ROFL_ERR("Error opening file")
                    ROFL_ASSERT(0)
                }

            SomUtils::MatD Tstart = SomUtils::MatD::Random(somSzD.d_ * somSzD.n_, 1);
            Tstart = startX.block(somSzD.d_ * somSzD.d_ * somSzD.n_, 0, somSzD.d_ * somSzD.n_, 1);
            Tstart.resize(somSzD.d_, somSzD.n_);

            ROFL_VAR1(startX)

            { // RSOM Procrustes execution scope
                rofl::ScopedTimer runProcrTimer("runRsomProcrustes");

                // method initialization and run
                SomProcrustes sp(somSzD, Tijs, edges);

                SomUtils::MatD Xgt = SomUtils::MatD::Random(somSzD.d_ * somSzD.d_ * somSzD.n_ + somSzD.d_ * somSzD.n_, 1);
                SomUtils::readCsvVecEigen(folderIn + "/Xgt.csv", Xgt);
                std::vector<SomUtils::MatD> Rgt(somSzD.n_, SomUtils::MatD::Identity(somSzD.d_, somSzD.d_));
                SomUtils::MatD RgtVec(SomUtils::MatD::Zero(somSzD.d_ * somSzD.d_ * somSzD.n_, 1));
                RgtVec = Xgt.block(0, 0, somSzD.d_ * somSzD.d_ * somSzD.n_, 1);
                SomUtils::MatD RgtHst = RgtVec.reshaped<Eigen::ColMajor>(somSzD.d_, somSzD.d_ * somSzD.n_);
                SomUtils::unStackH(RgtHst, Rgt, somSzD.d_);

                SomUtils::MatD Tgt = SomUtils::MatD::Random(somSzD.d_ * somSzD.n_, 1);
                Tgt = Xgt.block(somSzD.d_ * somSzD.d_ * somSzD.n_, 0, somSzD.d_ * somSzD.n_, 1);
                Tgt.resize(somSzD.d_, somSzD.n_);
                ROFL_VAR1(Tgt)

                // auto Tstart = Tgt;
                sp.setTstart(Tstart);

                sp.setTgt(Tgt); // TODO: add proper setter
                sp.setRgt(Rgt); // TODO: add proper setter

                sp.run();

                SomUtils::MatD Tout = sp.getTout();
                SomUtils::VecMatD Rout = sp.getRout();

                double costOut = sp.getCost();
                ROFL_VAR1(costOut)

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
                    startxofs << startX;
                }

                std::vector<double> rotErrsTestjd(numEdges, 1e+6), translErrsTestjd(numEdges, 1e+6);
                SomUtils::computeErrorsSingleRsom(edges,
                                                  Rout, Tout,
                                                  Rgt, Tgt,
                                                  rotErrsTestjd, translErrsTestjd);

                double execTimeIJ = runProcrTimer.elapsedTimeMs();
                ROFL_VAR3(entry, testjd, execTimeIJ)

                rotErrs[testjd].resize(numEdges, 1e+6);
                translErrs[testjd].resize(numEdges, 1e+6);

                rotErrs[testjd] = rotErrsTestjd;
                translErrs[testjd] = translErrsTestjd;

                execTimes[testjd] = execTimeIJ;

                for (int k = 0; k < numEdges; ++k)
                {
                    // output 2
                    // ROFL_VAR4(entry, j, rotErrs[i][j], translErrs[i][j]);
                    rotErrsOfstream << "testjd " + std::to_string(testjd) + " k " << std::to_string(k) << std::endl;
                    rotErrsOfstream << rotErrs[testjd][k] << std::endl;
                    translErrsOfstream << "testjd " + std::to_string(testjd) + " k " << std::to_string(k) << std::endl;
                    translErrsOfstream << translErrs[testjd][k] << std::endl;
                    ROFL_VAR2(rotErrs[testjd][k], translErrs[testjd][k])
                }
                execTimesOfstream << "j " + std::to_string(testjd) << std::endl;
                execTimesOfstream << execTimes[testjd] << std::endl;

                // Finding mean error of current instance-testjd pair
                rotMeanErr = SomUtils::stlVecDoublesMean(rotErrs[testjd]);
                translMeanErr = SomUtils::stlVecDoublesMean(translErrs[testjd]);
                rotErrsMeanOfstream << "j " + std::to_string(testjd) << std::endl;
                rotErrsMeanOfstream << rotMeanErr << std::endl;
                translErrsMeanOfstream << "j " + std::to_string(testjd) << std::endl;
                translErrsMeanOfstream << translMeanErr << std::endl;
                ROFL_VAR2(rotMeanErr, translMeanErr)

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
                ROFL_VAR2(rotErrs[j][k], translErrs[j][k])
            }
            ROFL_VAR2(rotMeanErr, translMeanErr)
        }

        // execTimesMeanOfstream << "i " + std::to_string(i) << std::endl;
        execTimesMeanOfstream << execTimeMean << std::endl;

        rotErrsAll.push_back(rotErrs);
        translErrsAll.push_back(translErrs);
        execTimesAll.push_back(execTimes);

        inst++; // current instance idx
    }

    rotErrsOfstream.close();
    translErrsOfstream.close();
    execTimesOfstream.close();
    rotErrsMeanOfstream.close();
    translErrsMeanOfstream.close();
    execTimesMeanOfstream.close();

    for (int i = 0; i < numInstances; ++i) // i already declared
    {
        for (int j = 0; j < numTestsPerInstance; ++j)
        {
            if (!rotErrsAll.empty())
            {
                for (int k = 0; k < rotErrsAll[i][j].size(); ++k)
                    ROFL_VAR5(i, j, k, rotErrsAll[i][j][k], translErrsAll[i][j][k]);
            }
            else
                continue;

            ROFL_VAR3(SomUtils::stlVecDoublesMean(rotErrsAll[i][j]), SomUtils::stlVecDoublesMean(translErrsAll[i][j]), SomUtils::stlVecDoublesMean(execTimesAll[i]));
        }
    }

    return 0;
}
