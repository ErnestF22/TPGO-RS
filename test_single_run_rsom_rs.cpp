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
    int d = 3;
    int n;
    SomUtils::readSingleIntCsv("../matlab/data/cpp_testdata/tdata_n10_mindeg3/n.csv", n);

    int numEdges;

    SomUtils::readSingleIntCsv("../matlab/data/cpp_testdata/tdata_n10_mindeg3/e.csv", numEdges);

    ROFL_VAR2(n, numEdges); // n = 10, numEdges = 48

    SomUtils::SomSize somSzD(d, d, n);
    SomUtils::MatD Tijs(d, numEdges);
    Eigen::MatrixXi edges(numEdges, 2);

    SomUtils::readMatlabCsvTijs("../matlab/data/cpp_testdata/tdata_n10_mindeg3/tijs.csv", Tijs, d, numEdges);
    SomUtils::readMatlabCsvEdges("../matlab/data/cpp_testdata/tdata_n10_mindeg3/edges.csv", edges);

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
    SomUtils::readCsvInitguess("../matlab/data/cpp_testdata/tdata_n10_mindeg3/X_gt.csv", xGt);
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

    bool readInitguess = true;
    if (readInitguess)
        SomUtils::readCsvInitguess("../data/X_initguess_test_single_run_rsom_rs.csv", startX);

    // RUN RSOM RS
    int srcNodeId = 0;
    ROPTLIB::runRsomRS(Prob, startX, srcNodeId); //note: startX is needed (even if random) in ROPTLIB; 
    //ROPTLIB namespace is used even if runRsomRS() is not in SampleSomProblem class, nor in "original" ROPTLIB


    
}