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

void removeQuasiZeros(Eigen::MatrixXd &mat, double thr = 1e-5)
{
    for (int i = 0; i < mat.rows(); ++i)
    {
        for (int j = 0; j < mat.cols(); ++j)
        {
            if (std::abs(mat(i, j)) < thr)
            {
                mat(i, j) = 0;
            }
        }
    }
}

SomUtils::MatD compareWithNans(const SomUtils::MatD &vec1, const SomUtils::MatD &vec2)
{

    const SomUtils::MatD zero_or_nan_1 = vec1 - vec1;
    const SomUtils::MatD zero_or_nan_2 = vec2 - vec2;

    const SomUtils::MatD compare = (vec1.array() < vec2.array()).matrix().cast<double>();

    return compare + zero_or_nan_1 + zero_or_nan_2;
}

bool compareWithNansEigenvectors(const SomUtils::MatD &vec, double &scalar, double tol = 1e-9)
{
    std::vector<double> vd;
    for (int i = 0; i < vec.rows(); ++i)
    {
        auto vecI = vec(i, 0);
        if (!std::isnan(vecI))
        {
            vd.push_back(vecI);
        }
    }

    const auto [minVd, maxVd] = std::minmax_element(std::begin(vd), std::end(vd));

    if (vd.size() > 0)
    {
        scalar = vd[0];
    }
    else
    {
        scalar = std::nan("nan");
    }

    ROFL_VAR2(*minVd, *maxVd)
    return (fabs(*maxVd - *minVd) < tol);
}

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
    SomUtils::SomSize somSzNext(d + 1, d, n);
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

    // SomUtils::VecMatD uRnext(somSzD.n_, SomUtils::MatD::Zero(somSzNext.p_, somSzNext.d_));
    // SomUtils::MatD uTnext(SomUtils::MatD::Zero(somSzNext.p_, somSzNext.n_));
    ProbNext.getRotations(Xnext, xRnext);
    ProbNext.getTranslations(Xnext, xTnext);
    // ProbNext.getRotations(Unext, uRnext);
    // ProbNext.getTranslations(Unext, uTnext);

    // problem_struct_next = problem_struct;
    // problem_struct_next.sz(1) = p;
    // Hgp = hess_genproc(Xnext,Unext,problem_struct_next);
    // disp("Hgp")
    // disp(Hgp)

    // SomUtils::VecMatD rhr(somSzD.n_, SomUtils::MatD::Zero(somSzNext.p_, somSzNext.d_));
    // SomUtils::MatD rht(SomUtils::MatD::Zero(somSzNext.p_, somSzNext.n_));
    // ProbNext.hessGenprocEigen(xRnext, uRnext, xTnext, uTnext, rhr, rht);
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
    std::ofstream hessianFile("hmat.csv");
    hessianFile << Hmat.reshaped<Eigen::ColMajor>(Hmat.cols() * Hmat.rows(), 1);
    hessianFile.close();

    if (SomUtils::isEqualFloats(Hmat - Hmat.transpose(), SomUtils::MatD::Zero(vecsz, vecsz)))
    {
        // Checking whether anti-symmetric part is zero
        ROFL_ERR("Hessian NOT symmetric")
        ROFL_ASSERT(0)
    }

    // Eigen::SelfAdjointEigenSolver<SomUtils::MatD> es;
    // es.compute(Hmat);
    // auto eigvals = es.eigenvalues(); // TODO: check how to use sparse matrix methods -> maybe SPECTRA library?
    // auto eigvecs = es.eigenvectors();
    // ROFL_VAR5(vecsz, eigvals.rows(), eigvals.cols(), eigvecs.rows(), eigvecs.cols());
    // ROFL_VAR1(Hmat)

    SomUtils::VecMatD xR(somSzD.n_, SomUtils::MatD::Zero(somSzD.d_, somSzD.d_));
    SomUtils::MatD xT(SomUtils::MatD::Zero(somSzD.d_, somSzD.n_));
    xT = xTnext.block(0, 0, somSzD.d_, somSzD.n_);
    for (int i = 0; i < somSzD.n_; ++i)
    {
        xR[i] = xRnext[i].block(0, 0, somSzD.d_, somSzD.d_);
    }
    SomUtils::VecMatD vRout(somSzD.n_, SomUtils::MatD::Zero(somSzNext.p_, somSzD.d_));
    SomUtils::MatD vTout(SomUtils::MatD::Zero(somSzNext.p_, somSzD.n_));
    double lambdaOut;
    Prob.rsomEscapeHessianGenprocEigen(xR, xT, Y0, lambdaOut, vRout, vTout);

    // vmin = V(:,x);
    // disp("vmin")
    // disp(vmin)

    ROFL_VAR1(lambdaOut)

    // [~, lambda_pim_out, v_pim_out] = rsom_pim_hessian_genproc(X, problem_struct);
    SomUtils::VecMatD vPimRout(somSzD.n_, SomUtils::MatD::Zero(somSzNext.p_, somSzD.d_));
    SomUtils::MatD vPimTout(SomUtils::MatD::Zero(somSzNext.p_, somSzD.n_));
    double lambdaPimOut;
    ROPTLIB::Vector Y0pim;
    Prob.rsomPimHessianGenprocEigen(1e-4, xR, xT, Y0pim, lambdaPimOut, vPimRout, vPimTout);

    ROFL_VAR1(lambdaPimOut)

    // disp("lambda_pim_out")
    // disp(lambda_pim_out)
    // disp("vectorizeXrt(v_pim_out)")
    // disp(vectorizeXrt(v_pim_out))
    // disp("vmin, vectorizeXrt(v_pim_out)")
    // disp([vmin, vectorizeXrt(v_pim_out)])

    SomUtils::MatD XvecPimOut = SomUtils::MatD::Zero(vecsz, 1);
    ProbNext.vectorizeRT(vPimRout, vPimTout, XvecPimOut);

    SomUtils::MatD XvecOut = SomUtils::MatD::Zero(vecsz, 1);
    ProbNext.vectorizeRT(vRout, vTout, XvecOut);

    // check whether vmin, v_pim_out are proportional

    double tol = 1e-9; // Divide tolerance
    // A = remove_quasi_zeros(vmin, tol);
    // B = remove_quasi_zeros(vectorizeXrt(v_pim_out), tol);
    // C = A./B; % Divide element-wise
    removeQuasiZeros(XvecPimOut, tol);
    removeQuasiZeros(XvecOut, tol);

    SomUtils::MatD tmp = SomUtils::MatD::Zero(vecsz, 2);
    tmp.col(0) = XvecPimOut;
    tmp.col(1) = XvecOut;
    ROFL_VAR1(tmp)
    auto A = XvecPimOut;
    auto B = XvecOut;
    auto C = A.cwiseQuotient(B); // Divide element-wise

    // removing NaNs (0-divisions)
    // id1 = isnan(C);
    // C(id1) = [];

    // % Check if the difference between largest and smallest values are within the tolerance
    // check = abs(max(C) - min(C)) < tol;
    double scalar = std::nan("nan");
    bool check = compareWithNansEigenvectors(C, scalar, tol);

    // If check OK, get the scalar multiple
    if (check)
    {
        ROFL_ERR("Eigenvectors check OK!")
    }
    else
    {
        ROFL_ERR("Eigenvectors check FAIL!")
    }

    // disp("scalar")
    // disp(scalar)
    ROFL_VAR1(scalar)

    // checking also eigenvalues
    ROFL_VAR2(lambdaOut, lambdaPimOut)

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