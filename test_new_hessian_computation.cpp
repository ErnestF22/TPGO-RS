#include <iostream>
#include <fstream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <filesystem>
#include <set>
#include <boost/lexical_cast.hpp>

#include "thirdparty/roptlib/problems_SampleSom.h"

#include <rofl/common/param_map.h>

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
        std::string("../matlab/data/cpp_testdata_noisy/tdata_n5_mindeg3_sigma00/"));

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

    std::string tijsPathStr = "../data/Tijs_tnhc.csv";
    if (!SomUtils::readMatlabCsvTijs(tijsPathStr, Tijs, d, numEdges))
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
    integer numoftypes = 2; // 2 i.e. (3D) Stiefel + Euclidean
    integer numofmani1 = n; // num of Stiefel manifolds
    integer numofmani2 = 1;
    ROPTLIB::Stiefel mani1(d, d);
    mani1.ChooseParamsSet2();
    ROPTLIB::Euclidean mani2(d, n);
    ROPTLIB::ProductManifold ProdMani(numoftypes, &mani1, numofmani1, &mani2, numofmani2);
    ROPTLIB::SampleSomProblem Prob(somSzD, Tijs, edges);
    SomUtils::SomSize somSzNext(d+1, d, n);
    ROPTLIB::SampleSomProblem ProbNext(somSzNext, Tijs, edges);


    std::string XnextPathStr = "../data/Xnext_vec.csv";
    std::string UnextPathStr = "../data/Unext_vec.csv";
    SomUtils::MatD Xnext(somSzD.n_ * somSzD.d_ * (somSzD.d_ + 1) + somSzD.n_ * (somSzD.d_ + 1), 1);
    if (!SomUtils::readCsvVecEigen(XnextPathStr, Xnext))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }

    // ROFL_VAR1(Xnext.transpose())

    SomUtils::MatD Unext(somSzD.n_ * somSzD.d_ * (somSzD.d_ + 1) + somSzD.n_ * (somSzD.d_ + 1), 1);
    if (!SomUtils::readCsvVecEigen(UnextPathStr, Unext))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }

    // ROFL_VAR1(Unext.transpose())


    SomUtils::VecMatD xRnext(somSzD.n_, SomUtils::MatD::Zero(somSzNext.p_, somSzNext.d_));
    SomUtils::MatD xTnext(SomUtils::MatD::Zero(somSzNext.p_, somSzNext.n_));
    ROPTLIB::Vector Y0;
    double lambdaPimOut;
    SomUtils::VecMatD vPimRout;
    SomUtils::MatD vPimTout;

    SomUtils::VecMatD uRnext(somSzD.n_, SomUtils::MatD::Zero(somSzNext.p_, somSzNext.d_));
    SomUtils::MatD uTnext(SomUtils::MatD::Zero(somSzNext.p_, somSzNext.n_));
    ProbNext.getRotations(Xnext, xRnext);
    ProbNext.getTranslations(Xnext, xTnext);
    ProbNext.getRotations(Unext, uRnext);
    ProbNext.getTranslations(Unext, uTnext);

    // problem_struct_next = problem_struct;
    // problem_struct_next.sz(1) = p;
    // Hgp = hess_genproc(Xnext,Unext,problem_struct_next);
    // disp("Hgp")
    // disp(Hgp)

    SomUtils::VecMatD rhr(somSzD.n_, SomUtils::MatD::Zero(somSzNext.p_, somSzNext.d_));
    SomUtils::MatD rht(SomUtils::MatD::Zero(somSzNext.p_, somSzNext.n_));
    ProbNext.hessGenprocEigen(xRnext, uRnext, xTnext, uTnext, rhr, rht);
    // for (auto i = 0; i < somSzD.n_; i++)
    // {
    //     ROFL_VAR1(rhr[i]);
    // }
    // ROFL_VAR1(rht);

    // Htest = Hmat * vectorizeXrt(Unext);
    // disp("[Hgp, Htest]")
    // disp([vectorizeXrt(Hgp), Htest])
    // disp("min(abs(vectorizeXrt(Hgp) - Htest), [], ""all"")")
    // disp(min(abs(vectorizeXrt(Hgp) - Htest), [], "all"))


    // [V, lambdas] = eig(Hmat);
    // % disp("V")
    // % disp(V)
    // % disp("lambdas")
    // % disp(lambdas)

    int vecsz = Xnext.rows();
    SomUtils::MatD Hmat(SomUtils::MatD::Zero(vecsz, vecsz));

    ProbNext.makeHmat(Xnext, somSzNext, Hmat);

    Eigen::EigenSolver<SomUtils::MatD> es;
    es.compute(Hmat);
    auto eigvals = es.eigenvalues(); // TODO: check how to use sparse matrix methods -> maybe SPECTRA library?
    auto eigvecs = es.eigenvectors();

    ROFL_VAR5(vecsz, eigvals.rows(), eigvals.cols(), eigvecs.rows(), eigvecs.cols());

    lambdaPimOut = eigvals.real().minCoeff();
    ROFL_VAR1(lambdaPimOut)

    // ROFL_VAR1(Hmat)

    // auto minEigvalIdx = eigvals.diagonal().argmin();
    // auto minEigvec = eigvecs.col(minEigvalIdx);
    // ProbNext.getRotations(minEigvec, vPimRout);
    // ProbNext.getTranslations(minEigvec, vPimTout);


    // Prob.rsomEscapeHessianGenprocEigen(xRnext, xTnext, Y0, lambdaPimOut, vPimRout, vPimTout);


    // max_imag_part = max(abs(imag(lambdas)), [], "all");

    // disp("max_imag_part")
    // disp(max_imag_part)

    // if (max_imag_part > 1e-6)
    //     error("Hessian NOT symmetric")
    // end

    // lambdas = real(lambdas);
    // % V = real(V);

    // lambda_min = min(diag(lambdas)); %!!
    // disp("lambda_min")
    // disp(lambda_min)

    // [x,y]=find(lambdas==lambda_min); %Note: x, y should be equal (lambdas matrix is diagonal)
    // disp("x")
    // disp(x)
    // disp("y")
    // disp(y)

    // vmin = V(:,x);
    // disp("vmin")
    // disp(vmin)



    // [~, lambda_pim_out, v_pim_out] = rsom_pim_hessian_genproc(X, problem_struct);

    // disp("lambda_pim_out")
    // disp(lambda_pim_out)
    // disp("vectorizeXrt(v_pim_out)")
    // disp(vectorizeXrt(v_pim_out))

    // disp("vmin, vectorizeXrt(v_pim_out)")
    // disp([vmin, vectorizeXrt(v_pim_out)])

    // % check whether vmin, v_pim_out are proportional

    // tol = 1e-10; % Divide tolerance
    // A = remove_quasi_zeros(vmin, 10 * tol);
    // B = remove_quasi_zeros(vectorizeXrt(v_pim_out), 10 * tol);
    // C = A./B; % Divide element-wise

    // %removing NaNs (0-divisions)
    // id1 = isnan(C);
    // C(id1) = [];

    // % Check if the difference between largest and smallest values are within the
    // % tolerance
    // check = abs(max(C) - min(C)) < tol;

    // % If yes, get the scalar multiple
    // if check
    //     scalar = C(1);
    //     disp("Eigenvectors check OK!")
    // else % If not, set to NaN
    //     scalar = NaN;
    // end

    // disp("scalar")
    // disp(scalar)

    // %checking also eigenvalue
    // disp("[lambda_min, lambda_pim_out]")
    // disp([lambda_min, lambda_pim_out])

    // lambdas_check = abs(lambda_min - lambda_pim_out) < 1e-3;
    // disp("abs(lambda_min - lambda_pim_out) < 1e-3")
    // disp(lambdas_check)

    // if lambdas_check
    //     disp("Eigenvalues check OK!")
    // else
    //     disp("Eigenvalues check FAIL!")
    // end

    return 0;
}