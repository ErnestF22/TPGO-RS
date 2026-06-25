#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <rofl/common/param_map.h>

#include "som_utils.h"
#include "thirdparty/roptlib/problems_Lsom.h"

void readVectorFromFile(const std::string &filename, SomUtils::MatD &vec);

void probRgradVectorized(ROPTLIB::LsomProblem &Prob, const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas, std::string outFilename);

void probRhessVectorized(ROPTLIB::LsomProblem &Prob, const SomUtils::VecMatD &R, const SomUtils::VecMatD &uR, const SomUtils::MatD &T, const SomUtils::MatD &uT, const SomUtils::MatD &Lambdas, const SomUtils::MatD &uLambdas, std::string outFilename);

int main(int argc, char **argv)
{
    std::string filenameCfg;
    std::string folderIn;

    int d = 3;
    int nrs = 3;
    int n = 5;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.getParam<std::string>("folderIn", folderIn, std::string("/home/rimlab/workspace/matlab_ws/som/matlab/2026_06_25_12_13_43/"));

    params.read(filenameCfg);
    params.read(argc, argv);

    int numEdges = 0;
    if (!SomUtils::readSingleIntCsv(folderIn + "/num_edges.csv", numEdges))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }
    ROFL_VAR1(numEdges)

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
    ROPTLIB::ProductManifold ProdManiLsom(numoftypes,
                                          &mani1, numofmani1, &mani2, numofmani2, &mani3, numofmani3);

    ROPTLIB::LsomProblem Prob(somSzD, Tijs, edges);

    double a = 0.0;
    if (!SomUtils::readSingleDoubleCsv(folderIn + "/a.csv", a))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }
    ROFL_VAR1(a)
    Prob.setLogScaleCompensationParam(a);

    double rho = 0.0;
    if (!SomUtils::readSingleDoubleCsv(folderIn + "/rho.csv", rho))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }
    ROFL_VAR1(rho)
    Prob.setRho(rho);

    double mu = 0.0;
    if (!SomUtils::readSingleDoubleCsv(folderIn + "/mu.csv", mu))
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }
    ROFL_VAR1(mu)

    Prob.setMuAdmm(mu);

    SomUtils::MatD y(numEdges, 1);
    readVectorFromFile(folderIn + "/y.csv", y);
    ROFL_VAR1(y)

    Prob.setYAdmm(y);

    SomUtils::MatD z(numEdges, 1);
    readVectorFromFile(folderIn + "/z.csv", z);
    ROFL_VAR1(z)

    Prob.setZAdmm(z);

    // Compute gradient and Hessian at point x, tangent point u

    SomUtils::MatD x(d * d * n + d * n + numEdges, 1);
    readVectorFromFile(folderIn + "/X_vec.csv", x);
    ROFL_VAR1(x)

    SomUtils::MatD u(d * d * n + d * n + numEdges, 1);
    readVectorFromFile(folderIn + "/X_tg_vec.csv", u);
    ROFL_VAR1(u)

    SomUtils::VecMatD R(n, SomUtils::MatD::Zero(d, d));
    Prob.getRotations(x, R);
    SomUtils::MatD T(d, n);
    Prob.getTranslations(x, T);
    SomUtils::MatD Lambdas(numEdges, 1);
    Prob.getScales(x, Lambdas);

    probRgradVectorized(Prob, R, T, Lambdas, folderIn + "/grad_X_cpp.csv");

    SomUtils::VecMatD uR(n, SomUtils::MatD::Zero(d, d));
    Prob.getRotations(u, uR);
    SomUtils::MatD uT(d, n);
    Prob.getTranslations(u, uT);
    SomUtils::MatD uLambdas(numEdges, 1);
    Prob.getScales(u, uLambdas);

    // ROFL_VAR1(Prob.hessGenprocEigen(x, u));

    probRhessVectorized(Prob, R, uR, T, uT, Lambdas, uLambdas, folderIn + "/hess_X_U_cpp.csv");

    return 0;
}

void probRgradVectorized(ROPTLIB::LsomProblem &Prob, const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas, std::string outFilename)
{
    Eigen::IOFormat OnePerLine(Eigen::StreamPrecision, Eigen::DontAlignCols, "\n", "\n", "", "", "", "");

    // Implement the vectorized version of the gradient computation

    int d = Prob.sz_.d_;
    int n = Prob.sz_.n_;
    int numEdges = Prob.numEdges_;

    SomUtils::VecMatD rgR(n, SomUtils::MatD::Zero(d, d));
    Prob.rgradR(R, T, Lambdas, rgR);
    SomUtils::MatD rgT(d, n);
    Prob.rgradT(R, T, Lambdas, rgT);
    SomUtils::MatD rgLambdas(numEdges, 1);
    Prob.rgradLambdas(R, T, Lambdas, rgLambdas);

    // Write the results to files
    std::ofstream gradOfs(outFilename);
    if (!gradOfs)
    {
        ROFL_ERR("Error opening file for writing gradient")
        ROFL_ASSERT(0)
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < d; ++j)
        {
            for (int k = 0; k < d; ++k)
            {
                gradOfs << rgR[i](k, j);
                gradOfs << std::endl;
            }
        }
    }
    gradOfs << rgT.reshaped(d * n, 1).format(OnePerLine) << std::endl;
    gradOfs << rgLambdas.format(OnePerLine) << std::endl;
}

void probRhessVectorized(ROPTLIB::LsomProblem &Prob, const SomUtils::VecMatD &R, const SomUtils::VecMatD &uR, const SomUtils::MatD &T, const SomUtils::MatD &uT, const SomUtils::MatD &Lambdas, const SomUtils::MatD &uLambdas, std::string outFilename)
{
    Eigen::IOFormat OnePerLine(Eigen::StreamPrecision, Eigen::DontAlignCols, "\n", "\n", "", "", "", "");

    
    // Implement the vectorized version of the Hessian computation

    int d = Prob.sz_.d_;
    int n = Prob.sz_.n_;
    int numEdges = Prob.numEdges_;

    SomUtils::VecMatD rhR(n, SomUtils::MatD::Zero(d, d));
    SomUtils::MatD rhT(d, n);
    SomUtils::MatD rhLambdas(numEdges, 1);

    Prob.hessGenprocEigen(R, uR, T, uT, Lambdas, uLambdas, rhR, rhT, rhLambdas);

    // Write the results to files
    std::ofstream hessOfs(outFilename);
    if (!hessOfs)
    {
        ROFL_ERR("Error opening file for writing Hessian")
        ROFL_ASSERT(0)
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < d; ++j)
        {
            for (int k = 0; k < d; ++k)
            {
                hessOfs << rhR[i](k, j);
                hessOfs << std::endl;
            }
        }
    }
    hessOfs << rhT.reshaped(d * n, 1).format(OnePerLine) << std::endl;
    hessOfs << rhLambdas.format(OnePerLine) << std::endl;
}

void readVectorFromFile(const std::string &filename, SomUtils::MatD &vec)
{

    std::ifstream ifs(filename);
    if (!ifs)
    {
        ROFL_ERR("Error opening file")
        ROFL_ASSERT(0)
    }

    std::vector<double> values;
    std::string line;
    while (std::getline(ifs, line))
    {
        // Allow comma or whitespace separated numbers; skip empty lines
        if (line.find_first_not_of(" \t\r\n") == std::string::npos)
            continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream ss(line);
        double v;
        while (ss >> v)
            values.push_back(v);
    }

    for (size_t i = 0; i < values.size(); ++i)
        vec(static_cast<Eigen::Index>(i), 0) = values[i];

    if (static_cast<int>(values.size()) != vec.rows())
    {
        ROFL_ERR("vector size does not match numEdges")
        ROFL_ASSERT(0)
    }

    // ROFL_VAR1(vec)
}
