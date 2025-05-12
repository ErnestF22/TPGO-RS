#include <iostream>
#include <fstream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <filesystem>
#include <boost/lexical_cast.hpp>

#include "thirdparty/roptlib/problems_SampleSom.h"

int main(int argc, char **argv)
{
    int p = 4;
    int d = 3;
    int n = 10;

    double alpha = 1.0;

    int numEdges = 48; // should be unused

    SomUtils::SomSize somSz(p, d, n);
    SomUtils::MatD Tijs(d, numEdges);
    Eigen::MatrixXi edges(numEdges, 2);

    SomUtils::readMatlabCsvTijs("../matlab/data/cpp_testdata/tdata_n10_mindeg3/tijs.csv", Tijs, d, numEdges);
    SomUtils::readMatlabCsvEdges("../matlab/data/cpp_testdata/tdata_n10_mindeg3/edges.csv", edges);

    integer numoftypes = 2; // 2 i.e. (3D) Stiefel + Euclidean
    integer numofmani1 = n; // num of Stiefel manifolds
    integer numofmani2 = 1;
    ROPTLIB::Stiefel mani1(p, d);
    mani1.ChooseParamsSet2();
    ROPTLIB::Euclidean mani2(p, n);
    ROPTLIB::ProductManifold ProdMani(numoftypes, &mani1, numofmani1, &mani2, numofmani2);
    ROPTLIB::SampleSomProblem Prob(somSz, Tijs, edges);

    // x
    SomUtils::MatD xVecEig(SomUtils::MatD::Zero(p * d * n + p * n, 1));
    SomUtils::VecMatD REig(n, SomUtils::MatD::Zero(p, d));
    SomUtils::MatD TEig(SomUtils::MatD::Zero(p, n));
    ROPTLIB::Vector x = ProdMani.RandominManifold();
    SomUtils::readCsvInitguess("../matlab/data/retraction_debug/matlab_x.csv", x);
    Prob.RoptToEig(x, xVecEig);
    ROFL_VAR2(xVecEig.rows(), xVecEig.cols())
    Prob.getRotations(xVecEig, REig);
    Prob.getTranslations(xVecEig, TEig);

    for (int i = 0; i < n; ++i)
        ROFL_VAR2(i, REig[i])

    ROFL_VAR1(TEig)

    // d
    SomUtils::MatD dVecEig(SomUtils::MatD::Zero(p * d * n + p * n, 1));
    SomUtils::VecMatD dREig(n, SomUtils::MatD::Zero(p, d));
    SomUtils::MatD dTEig(SomUtils::MatD::Zero(p, n));
    ROPTLIB::Vector u = ProdMani.RandominManifold();
    SomUtils::readCsvInitguess("../matlab/data/retraction_debug/matlab_d.csv", u);
    Prob.RoptToEig(u, dVecEig);
    ROFL_VAR2(dVecEig.rows(), dVecEig.cols())
    Prob.getRotations(dVecEig, dREig);
    Prob.getTranslations(dVecEig, dTEig);

    for (int i = 0; i < n; ++i)
        ROFL_VAR2(i, dREig[i])

    ROFL_VAR1(dTEig)
    ROFL_VAR1("\n")

    // newx
    SomUtils::MatD newxMlabVecEig(SomUtils::MatD::Zero(p * d * n + p * n, 1));
    SomUtils::VecMatD newxMlabREig(n, SomUtils::MatD::Zero(p, d));
    SomUtils::MatD newxMlabTEig(SomUtils::MatD::Zero(p, n));
    ROPTLIB::Vector newxMlab = ProdMani.RandominManifold();
    SomUtils::readCsvInitguess("../matlab/data/retraction_debug/matlab_newx.csv", newxMlab);
    Prob.RoptToEig(newxMlab, newxMlabVecEig);
    ROFL_VAR2(newxMlabVecEig.rows(), newxMlabVecEig.cols())
    Prob.getRotations(newxMlabVecEig, newxMlabREig);
    Prob.getTranslations(newxMlabVecEig, newxMlabTEig);

    // test retraction
    SomUtils::VecMatD newxCppR(n, SomUtils::MatD::Zero(p, d));
    SomUtils::MatD newxCppT(SomUtils::MatD::Zero(p, n));
    Prob.stiefelRetractionQR(REig, dREig, newxCppR, alpha);
    Prob.euclRetraction(TEig, dTEig, newxCppT, alpha);

    // R - QR
    for (int i = 0; i < n; ++i)
        ROFL_VAR3(i, newxCppR[i], newxMlabREig[i])

    // T
    ROFL_VAR2(newxCppT, newxMlabTEig)

    // R - Polar
    ROFL_VAR1("POLAR Stiefel retr")
    SomUtils::VecMatD newxRpolarMlab(n, SomUtils::MatD::Zero(p, d));
    SomUtils::MatD newxRpolarMlabHst(SomUtils::MatD::Zero(p, d * n));
    SomUtils::readCsvEigen("../matlab/data/retraction_debug/matlab_newxRpolar.csv", newxRpolarMlabHst);
    SomUtils::unStackH(newxRpolarMlabHst, newxRpolarMlab);
    Prob.stiefelRetractionPolar(REig, dREig, newxCppR, alpha);
    for (int i = 0; i < n; ++i)
        ROFL_VAR3(i, newxRpolarMlab[i], newxMlabREig[i])

    return 0;
}