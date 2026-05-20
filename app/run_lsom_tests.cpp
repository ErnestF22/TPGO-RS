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

#include "thirdparty/roptlib/problems_Lsom.h"

#include "som_utils.h"

namespace fs = std::filesystem;

void EigToRopt(const SomUtils::MatD &xEig, const ROPTLIB::LsomProblem &Prob, ROPTLIB::Vector *result);

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
    bool firstZadmmEqualToLambdas;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>(
        "in", folderIn,
        std::string("../matlab/data/ssom_testdata_noisy/small_sample/"));

    params.getParam<int>("d", d, 3);
    params.getParam<int>("numTestsPerInstance", numTestsPerInstance, 50);
    params.getParam<bool>("readStartingPtFromFile", readStartingPtFromFile, false);
    params.getParam<std::string>("resultsBasePath", resultsBasePath, "../results_lsom/");
    params.getParam<int>("srcNodeIdx", srcNodeIdx, 0);
    params.getParam<double>("rho", rho, 1000.0);
    params.getParam<bool>("firstZadmmEqualToLambdas", firstZadmmEqualToLambdas, false);

    std::cout << "Params:" << std::endl;
    params.write(std::cout);

    std::cout << "-------\n"
              << std::endl;

    /*************************End of ROFL params reading**************************/

    for (auto &entry : fs::directory_iterator(folderIn))
        sortedByName.insert(entry.path());

    std::vector<std::vector<std::vector<double>>> rotErrsAll, translErrsAll, lambdaErrsAll;
    std::vector<std::vector<double>> execTimesAll;
    std::vector<std::vector<double>> staircaseStepOutIdxAll;

    // declaring ofstreams
    std::ofstream rotErrsOfstream;
    std::ofstream translErrsOfstream;
    std::ofstream lambdaErrsOfstream;
    std::ofstream execTimesOfstream;
    std::ofstream staircaseStepOutIdxOfstream;
    std::ofstream rsSuccessOfstream;
    std::ofstream rotDetsOkOfstream;
    std::ofstream lambdasAcceptableOfstream;
    std::ofstream rotErrsMeanOfstream;
    std::ofstream translErrsMeanOfstream;
    std::ofstream lambdaErrsMeanOfstream;
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
            ROFL_VAR1(entry.string() + "/n.csv")
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
        SomUtils::MatD tijs(d, numEdges);
        Eigen::MatrixXi edges(numEdges, 2);

        if (!SomUtils::readMatlabCsvTijs(entry.string() + "/tijs.csv", tijs, d, numEdges))
        {
            ROFL_ERR("Error opening file")
            ROFL_ASSERT(0)
        }
        if (!SomUtils::readMatlabCsvEdges(entry.string() + "/edges.csv", edges))
        {
            ROFL_ERR("Error opening file")
            ROFL_ASSERT(0)
        }

        ROFL_VAR1(tijs)
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
        ROPTLIB::ProductManifold ProdManiLsom(numoftypes,
                                              &mani1, numofmani1, &mani2, numofmani2, &mani3, numofmani3);

        SomUtils::MatD tijsNois = tijs;
        for (int e = 0; e < numEdges; ++e)
        {
            tijsNois.col(e) = tijs.col(e) / tijs.col(e).norm();
        }
        ROPTLIB::LsomProblem Prob(somSzD, tijsNois, edges);

        // Read GT from csv
        ROPTLIB::Vector xGt = ProdManiLsom.RandominManifold();
        if (!SomUtils::readCsvInitguess(entry.string() + "/Xgt.csv", xGt))
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
        Prob.setRho(rho); // default 1000.0
        Prob.RoptToEig(xGt, XgtVecEig);
        Prob.getRotations(XgtVecEig, RgtEig);
        Prob.getTranslations(XgtVecEig, TgtEig);
        Prob.getScales(XgtVecEig, LambdasGtEig);
        Prob.setGt(RgtEig, TgtEig, LambdasGtEig);

        // Set the domain of the problem to be the product of Stiefel manifolds
        Prob.SetDomain(&ProdManiLsom);

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
        // R
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
        // T
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
        // Lambda
        std::string lambdaErrsFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_lambda_errors.txt";
        if (lambdaErrsOfstream.is_open())
        {
            lambdaErrsOfstream.close();
            lambdaErrsOfstream.clear(); // clear flags
        }
        lambdaErrsOfstream.open(lambdaErrsFilename);
        if (!lambdaErrsOfstream)
        {
            ROFL_ERR("Error opening output file")
            ROFL_ASSERT(0)
        }
        // exec times
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
        // staircase step idx out
        std::string staircaseStepOutIdxFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_last_rs_step.txt";
        if (staircaseStepOutIdxOfstream.is_open())
        {
            staircaseStepOutIdxOfstream.close();
            staircaseStepOutIdxOfstream.clear(); // clear flags
        }
        staircaseStepOutIdxOfstream.open(staircaseStepOutIdxFilename);
        if (!staircaseStepOutIdxOfstream)
        {
            ROFL_ERR("Error opening output file")
            ROFL_ASSERT(0)
        }
        // rs success
        std::string rsSuccessFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_rs_success.txt";
        if (rsSuccessOfstream.is_open())
        {
            rsSuccessOfstream.close();
            rsSuccessOfstream.clear(); // clear flags
        }
        rsSuccessOfstream.open(rsSuccessFilename);
        if (!rsSuccessOfstream)
        {
            ROFL_ERR("Error opening output file")
            ROFL_ASSERT(0)
        }
        // rot dets ok
        std::string rotDetsOkFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_rot_dets_ok.txt";
        if (rotDetsOkOfstream.is_open())
        {
            rotDetsOkOfstream.close();
            rotDetsOkOfstream.clear(); // clear flags
        }
        rotDetsOkOfstream.open(rotDetsOkFilename);
        if (!rotDetsOkOfstream)
        {
            ROFL_ERR("Error opening output file")
            ROFL_ASSERT(0)
        }
        // lambdas acceptable
        std::string lambdasAcceptableFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_lambdas_acceptable.txt";
        if (lambdasAcceptableOfstream.is_open())
        {
            lambdasAcceptableOfstream.close();
            lambdasAcceptableOfstream.clear(); // clear flags
        }
        lambdasAcceptableOfstream.open(lambdasAcceptableFilename);
        if (!lambdasAcceptableOfstream)
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
        std::string lambdaErrsMeanFilename = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_lambda_errors_mean.txt";
        if (lambdaErrsMeanOfstream.is_open())
        {
            lambdaErrsMeanOfstream.close();
            lambdaErrsMeanOfstream.clear(); // clear flags
        }
        lambdaErrsMeanOfstream.open(lambdaErrsMeanFilename);
        if (!lambdaErrsMeanOfstream)
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
        double rotMeanErr = 1e+6, translMeanErr = 1e+6, lambdaMeanErr = 1e+6, execTimeMean = 1e+6;
        std::vector<std::vector<double>> rotErrs(numTestsPerInstance), translErrs(numTestsPerInstance), lambdaErrs(numTestsPerInstance);
        std::vector<double> execTimes(numTestsPerInstance), staircaseStepOutIdx(numTestsPerInstance);
        std::vector<bool> rsSuccess(numTestsPerInstance), rotDetsOk(numTestsPerInstance), lambdasAcceptable(numTestsPerInstance);

        for (int testjd = 0; testjd < numTestsPerInstance; ++testjd)
        {
            ROFL_VAR3(entry, testjd, "start");

            // Generate startX (random)
            ROPTLIB::Vector startX = ProdManiLsom.RandominManifold();
            // startX.Initialization(numoftypes, &mani1, numofmani1, &mani2, numofmani2);

            if (readStartingPtFromFile)
                if (!SomUtils::readCsvInitguess(entry.string() + "/ssom_x_start.csv", startX))
                {
                    ROFL_ERR("Error opening file")
                    ROFL_ASSERT(0)
                }

            startX.Print("startX");

            SomUtils::MatD startXeig = SomUtils::MatD::Zero(d * d * n + d * n + numEdges, 1);
            Prob.RoptToEig(startX, startXeig);
            SomUtils::MatD scalesInitguess = SomUtils::MatD::Ones(numEdges, 1);
            if (Prob.reluScaleCompensation_)
                scalesInitguess = 5.0 * SomUtils::MatD::Ones(numEdges, 1);
            else
                scalesInitguess = 10.0 * SomUtils::MatD::Ones(numEdges, 1);
            startXeig.block(d * d * n + d * n, 0, numEdges, 1) = scalesInitguess;

            Prob.setFirstZadmmLambdas(firstZadmmEqualToLambdas);

            ROPTLIB::Vector startX2 = ProdManiLsom.RandominManifold();
            EigToRopt(startXeig, Prob, &startX2);

            startX2.Print("startX2");

            Prob.setPerformGlobalization(true);

            { // ssom RS execution scope

                rofl::ScopedTimer runSsomRStimer("runSsomRS");

                // RUN ssom RS
                SomUtils::VecMatD Rout(n, SomUtils::MatD::Identity(d, d));
                SomUtils::MatD Tout(SomUtils::MatD::Zero(d, n));
                SomUtils::MatD LambdasOut(SomUtils::MatD::Zero(numEdges, 1));
                int staircaseStepIdxOutIJ;
                bool rsSuccessIJ, rotDetsOkIJ, lambdasAcceptableIJ;
                /* Setting up Prob using setters */
                Prob.setUsePIM(true);           // same as default
                Prob.setPimMaxIterations(5000); // same as default

                rofl::ScopedTimer timer("ssomRS");

                double costOut = ROPTLIB::runLsom(Prob, startX2, srcNodeIdx, // !! startX2 used here!
                                                  Rout, Tout, LambdasOut,
                                                  staircaseStepIdxOutIJ,
                                                  rsSuccessIJ, rotDetsOkIJ, lambdasAcceptableIJ); // note: startX is needed (even if random) in ROPTLIB;

                // ROPTLIB namespace is used even if runRsomRS() is not in LsomProblem class, nor in "original" ROPTLIB
                ROFL_VAR1(costOut)
                // ROPTLIB namespace is used even if runssomRS() is not in LsomProblem class, nor in "original" ROPTLIB

                // costOut = 1.0f;
                if (!SomUtils::isEqualDoubles(costOut, 0.0f))
                {
                    auto fstr = resultsBasePath + folderAppendName + "_" + folderAppendNameStamped + "/" + folderAppendName + "_j" + std::to_string(testjd);
                    ROFL_VAR1(costOut)
                    std::ofstream tijsofs(fstr + "_tijs.txt");
                    tijsofs << tijs;
                    std::ofstream edgesofs(fstr + "_edges.txt");
                    edgesofs << edges;
                    std::ofstream lambdasOutOfs(fstr + "_lambdas_out.txt");
                    lambdasOutOfs << LambdasOut;
                    // SomUtils::MatD startXeig(SomUtils::MatD::Zero(d * d * n + d * n, 1));
                    // Prob.RoptToEig(startX, startXeig);
                    std::ofstream startxofs(fstr + "_startx.txt");
                    startxofs << startXeig; // startXeig is Eig version of startX2
                    tijsofs.close();
                    edgesofs.close();
                    lambdasOutOfs.close();
                    startxofs.close();
                }

                std::vector<double> rotErrsTestjd(numEdges, 1e+6), translErrsTestjd(numEdges, 1e+6), lambdaErrsTestjd(numEdges, 1e+6);
                SomUtils::computeErrorsSingleSsom(edges,
                                                  Rout, Tout, LambdasOut,
                                                  RgtEig, TgtEig, LambdasGtEig,
                                                  rotErrsTestjd, translErrsTestjd, lambdaErrsTestjd);

                double execTimeIJ = runSsomRStimer.elapsedTimeMs();
                ROFL_VAR3(entry, testjd, execTimeIJ)

                rotErrs[testjd].resize(numEdges, 1e+6);
                translErrs[testjd].resize(numEdges, 1e+6);
                lambdaErrs[testjd].resize(numEdges, 1e+6);

                rotErrs[testjd] = rotErrsTestjd;
                translErrs[testjd] = translErrsTestjd;
                lambdaErrs[testjd] = lambdaErrsTestjd;

                execTimes[testjd] = execTimeIJ;
                staircaseStepOutIdx[testjd] = staircaseStepIdxOutIJ;
                rsSuccess[testjd] = rsSuccessIJ;
                rotDetsOk[testjd] = rotDetsOkIJ;
                lambdasAcceptable[testjd] = lambdasAcceptableIJ;

                for (int k = 0; k < numEdges; ++k)
                {
                    // output 2
                    // ROFL_VAR4(entry, j, rotErrs[i][j], translErrs[i][j]);
                    rotErrsOfstream << "testjd " + std::to_string(testjd) + " k " << std::to_string(k) << std::endl;
                    rotErrsOfstream << rotErrs[testjd][k] << std::endl;
                    translErrsOfstream << "testjd " + std::to_string(testjd) + " k " << std::to_string(k) << std::endl;
                    translErrsOfstream << translErrs[testjd][k] << std::endl;
                    lambdaErrsOfstream << "testjd " + std::to_string(testjd) + " k " << std::to_string(k) << std::endl;
                    lambdaErrsOfstream << lambdaErrs[testjd][k] << std::endl;
                    ROFL_VAR2(rotErrs[testjd][k], translErrs[testjd][k])
                }
                execTimesOfstream << "j " + std::to_string(testjd) << std::endl;
                execTimesOfstream << execTimes[testjd] << std::endl;
                staircaseStepOutIdxOfstream << "j " + std::to_string(testjd) << std::endl;
                staircaseStepOutIdxOfstream << staircaseStepOutIdx[testjd] << std::endl;
                rsSuccessOfstream << "j " + std::to_string(testjd) << std::endl;
                rsSuccessOfstream << rsSuccess[testjd] << std::endl;
                rotDetsOkOfstream << "j " + std::to_string(testjd) << std::endl;
                rotDetsOkOfstream << rotDetsOk[testjd] << std::endl;
                lambdasAcceptableOfstream << "j " + std::to_string(testjd) << std::endl;
                lambdasAcceptableOfstream << lambdasAcceptable[testjd] << std::endl;

                // Finding mean error of current instance-testjd pair
                rotMeanErr = SomUtils::stlVecDoublesMean(rotErrs[testjd]);
                translMeanErr = SomUtils::stlVecDoublesMean(translErrs[testjd]);
                lambdaMeanErr = SomUtils::stlVecDoublesMean(lambdaErrs[testjd]);
                rotErrsMeanOfstream << "j " + std::to_string(testjd) << std::endl;
                rotErrsMeanOfstream << rotMeanErr << std::endl;
                translErrsMeanOfstream << "j " + std::to_string(testjd) << std::endl;
                translErrsMeanOfstream << translMeanErr << std::endl;
                lambdaErrsMeanOfstream << "j " + std::to_string(testjd) << std::endl;
                lambdaErrsMeanOfstream << lambdaMeanErr << std::endl;
                ROFL_VAR3(rotMeanErr, translMeanErr, lambdaMeanErr)

                execTimeMean = SomUtils::stlVecDoublesMean(execTimes);
            } // end of ssom RS execution scope
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
                ROFL_VAR3(rotErrs[j][k], translErrs[j][k], lambdaErrs[j][k])
            }
            ROFL_VAR3(rotMeanErr, translMeanErr, lambdaMeanErr)
        }

        // execTimesMeanOfstream << "i " + std::to_string(i) << std::endl;
        execTimesMeanOfstream << execTimeMean << std::endl;

        rotErrsAll.push_back(rotErrs);
        translErrsAll.push_back(translErrs);
        lambdaErrsAll.push_back(lambdaErrs);
        execTimesAll.push_back(execTimes);

        inst++; // current instance idx
    }

    rotErrsOfstream.close();
    translErrsOfstream.close();
    lambdaErrsOfstream.close();
    execTimesOfstream.close();
    staircaseStepOutIdxOfstream.close();
    rsSuccessOfstream.close();
    rotDetsOkOfstream.close();
    lambdasAcceptableOfstream.close();
    rotErrsMeanOfstream.close();
    translErrsMeanOfstream.close();
    lambdaErrsMeanOfstream.close();
    execTimesMeanOfstream.close();

    for (int i = 0; i < numInstances; ++i) // i already declared
    {
        for (int j = 0; j < numTestsPerInstance; ++j)
        {
            if (!rotErrsAll.empty())
            {
                for (int k = 0; k < rotErrsAll[i][j].size(); ++k)
                    ROFL_VAR6(i, j, k, rotErrsAll[i][j][k], translErrsAll[i][j][k], lambdaErrsAll[i][j][k]);
            }
            else
                continue;

            ROFL_VAR4(SomUtils::stlVecDoublesMean(rotErrsAll[i][j]), SomUtils::stlVecDoublesMean(translErrsAll[i][j]), SomUtils::stlVecDoublesMean(lambdaErrsAll[i][j]), SomUtils::stlVecDoublesMean(execTimesAll[i]));
        }
    }

    return 0;
}

void EigToRopt(const SomUtils::MatD &xEig, const ROPTLIB::LsomProblem &Prob, ROPTLIB::Vector *result)
{
    int rotSz = Prob.getRotSz();
    int translSz = Prob.getTranslSz();

    auto rgR = SomUtils::VecMatD(Prob.sz_.n_, SomUtils::MatD::Zero(Prob.sz_.p_, Prob.sz_.d_));
    // ROFL_VAR3(xEig.rows(), R.size(), T.size());
    // ROFL_VAR3(R.size(), P.size(), rgR.size());
    Prob.getRotations(xEig, rgR);
    SomUtils::MatD rgT(SomUtils::MatD::Zero(Prob.sz_.p_, Prob.sz_.n_));
    Prob.getTranslations(xEig, rgT);
    SomUtils::MatD rgLambdas(SomUtils::MatD::Zero(Prob.numEdges_, 1));
    Prob.getScales(xEig, rgLambdas);

    auto sz = Prob.sz_;

    int gElemIdx = 0;
    // fill result with computed gradient values: R
    for (int i = 0; i < sz.n_; ++i)
    {
        // ROFL_VAR1(gElemIdx);
        // ROFL_VAR2("\n", rgR[gElemIdx]);
        // result->GetElement(gElemIdx).SetToIdentity(); // Ri
        // result->GetElement(gElemIdx).Print("Ri before assignment");

        ROPTLIB::Vector rgRiVec(sz.p_, sz.d_);
        // rgRiVec.Initialize();
        realdp *GroptlibWriteArray = rgRiVec.ObtainWriteEntireData();
        for (int j = 0; j < rotSz; ++j)
        {
            // ROFL_VAR2(i, j);
            // rgRiVec.Print("rgRiVec before assignment");

            // ROFL_VAR1(rgRiVec.GetElement(j, 0));

            GroptlibWriteArray[j] = rgR[i].reshaped(sz.d_ * sz.p_, 1)(j);

            // ROFL_VAR1("");
            // rgRiVec.Print("rgRiVec after assignment");
        }
        rgRiVec.CopyTo(result->GetElement(gElemIdx));
        // result->GetElement(gElemIdx).Print("Riem. grad Ri after assignment");
        gElemIdx++;
    }

    // fill result with computed gradient values: T

    ROPTLIB::Vector rgTiVec(sz.p_, sz.n_);
    realdp *GroptlibWriteArray = rgTiVec.ObtainWriteEntireData();
    for (int j = 0; j < sz.p_ * sz.n_; ++j)
    {
        // rgTiVec.Print("rgTiVec before assignment");

        // ROFL_VAR1(rgRiVec.GetElement(j, 0));

        GroptlibWriteArray[j] = rgT.reshaped(sz.n_ * sz.p_, 1)(j);

        // ROFL_VAR1("");
        // rgTiVec.Print("rgTiVec after assignment");
    }
    rgTiVec.CopyTo(result->GetElement(gElemIdx));
    gElemIdx++;
    // result->GetElement(gElemIdx).Print("grad Ti after assignment");

    // ROFL_VAR2("\n", rgT);

    // fill result with computed gradient values: Lambdas

    ROPTLIB::Vector rgLambdasIvec(Prob.numEdges_, 1);
    realdp *GroptlibWriteArray2 = rgLambdasIvec.ObtainWriteEntireData();
    for (int j = 0; j < Prob.numEdges_; ++j)
    {
        // rhTiVec.Print("rhTiVec before assignment");

        // ROFL_VAR1(rhRiVec.GetElement(j, 0));

        GroptlibWriteArray2[j] = rgLambdas(j); // TODO: reshaped() call can probably be removed

        // ROFL_VAR1("");
        // rhTiVec.Print("rhTiVec after assignment");
    }
    rgLambdasIvec.CopyTo(result->GetElement(gElemIdx));
}