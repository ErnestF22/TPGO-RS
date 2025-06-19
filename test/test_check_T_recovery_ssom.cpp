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

struct procrustesOutput
{
    SomUtils::MatD T;
    double b;
    SomUtils::MatD c;

    procrustesOutput(const SomUtils::MatD &T_, double b_, const SomUtils::MatD &c_)
        : T(T_), b(b_), c(c_)
    {
    }
    procrustesOutput() : T(SomUtils::MatD::Zero(1, 1)), b(0.0), c(SomUtils::MatD::Zero(1, 1)) {}

    procrustesOutput(const procrustesOutput &other)
        : T(other.T), b(other.b), c(other.c)
    {
    }

    ~procrustesOutput() = default;
};

procrustesOutput procrustesCodemeta(const SomUtils::MatD &X, const SomUtils::MatD &Y, SomUtils::MatD &T, bool doScaling = true, const std::string &doReflection = "best");

void devectorizeXeig(const SomUtils::MatD &xEig, SomUtils::VecMatD &R, SomUtils::MatD &T, SomUtils::MatD &Lambdas, int p, int d, int n, int numEdges);

int main(int argc, char **argv)
{
    std::string filenameCfg;
    std::string folderIn;
    std::set<fs::path> sortedByName;

    int d;
    // int numTestsPerInstance;
    // bool readStartingPtFromFile;s
    int srcNodeIdx;

    rofl::ParamMap params;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    params.read(filenameCfg);
    params.read(argc, argv);
    // Output mode (quat or aa)
    params.getParam<std::string>(
        "in", folderIn,
        std::string("../data/check_T_recovery_ssom/"));

    params.getParam<int>("d", d, 3);
    // params.getParam<int>("numTestsPerInstance", numTestsPerInstance, 30);
    // params.getParam<bool>("readStartingPtFromFile", readStartingPtFromFile, true);
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

    int p;
    if (!SomUtils::readSingleIntCsv(folderIn + "p.csv", p))
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
    // ROPTLIB::Vector startX = ProdManiSsom.RandominManifold();

    SomUtils::MatD XmanoptOut(SomUtils::MatD::Zero(p * d * n + p * n + numEdges, 1));
    SomUtils::VecMatD RmanoptOut(n, SomUtils::MatD::Zero(p, d));
    SomUtils::MatD TmanoptOut(SomUtils::MatD::Zero(p, n));
    SomUtils::MatD LambdasManoptOut(SomUtils::MatD::Zero(numEdges, 1));

    if (!SomUtils::readCsvVecEigen(folderIn + "X_manopt_out.csv", XmanoptOut))
    {
        // matlab/data/ssom_testdata_noisy/harder/tdata_n5_mindeg2_sigma00/ssom_x_start.csv
        ROFL_ERR("Error opening file")
        ROFL_VAR1(folderIn + "X_manopt_out.csv")
        ROFL_ASSERT(0)
    }
    ROFL_VAR1(XmanoptOut.transpose());

    devectorizeXeig(XmanoptOut, RmanoptOut, TmanoptOut, LambdasManoptOut, p, d, n, numEdges);
    ROFL_VAR1("Printing R, T manopt out")
    for (auto &m : RmanoptOut)
        ROFL_VAR1(m)
    ROFL_VAR1(TmanoptOut)
    ROFL_VAR1(LambdasManoptOut.transpose())

    // end of test data reading

    SomUtils::MatD Tedges(SomUtils::MatD::Zero(p, numEdges));
    Prob.makeTedges(TmanoptOut, Tedges);
    ROFL_VAR1(Tedges)

    // RT_stacked_high_deg = [matStackH(R_manopt_out(:,:,nodes_high_deg)), T_edges];
    Eigen::ArrayXi nodeDegrees(Eigen::ArrayXi::Zero(n, 1));
    Prob.computeNodeDegrees(nodeDegrees);
    ROFL_VAR1(nodeDegrees)

    int lowDeg = 2; // low degree threshold
    Eigen::ArrayXi nodesHighDeg(Eigen::ArrayXi::Zero(n));
    Eigen::ArrayXi nodesLowDeg(Eigen::ArrayXi::Zero(n));

    for (int i = 0; i < n; ++i)
    {
        if (nodeDegrees(i, 0) > lowDeg)
            nodesHighDeg(i, 0) = 1;
        else
            nodesLowDeg(i, 0) = 1;
    }

    auto numNodesHighDeg = nodesHighDeg.sum();
    auto numNodesLowDeg = n - numNodesHighDeg;
    ;

    // Qalign = align3d(RT_stacked_high_deg);
    SomUtils::MatD RTstackedHighDeg(SomUtils::MatD::Zero(p, d * numNodesHighDeg + numEdges));
    SomUtils::MatD RstackedHighDeg(SomUtils::MatD::Zero(p, d * numNodesHighDeg));
    SomUtils::VecMatD RmanoptOutHighDeg;
    for (int i = 0; i < n; ++i)
    {
        if (nodesHighDeg(i, 0) != 0)
            RmanoptOutHighDeg.push_back(RmanoptOut[i]);
    }
    ROFL_VAR1("hstack call from here");
    SomUtils::hstack(RmanoptOutHighDeg, RstackedHighDeg);
    RTstackedHighDeg.block(0, 0, p, d * numNodesHighDeg) = RstackedHighDeg;
    RTstackedHighDeg.block(0, d * numNodesHighDeg, p, numEdges) = Tedges;
    SomUtils::MatD Qalign(SomUtils::MatD::Zero(p, p));
    Prob.align3d(RTstackedHighDeg, Qalign);
    ROFL_VAR1(Qalign)

    ROFL_VAR2(Qalign * RTstackedHighDeg, (Qalign * RTstackedHighDeg).block(d, 0, p - d, d * numNodesHighDeg + numEdges).cwiseAbs().maxCoeff());

    // Tij_2deg_recovery = [];
    // Tij_tilde_2deg_recovery = [];
    SomUtils::VecMatD Tij2degRecovery(numNodesLowDeg, SomUtils::MatD::Zero(p, d));
    SomUtils::VecMatD TijTilde2degRecovery(numNodesLowDeg, SomUtils::MatD::Zero(p, d));
    // for node_id = 1:N
    //     if problem_data.node_degrees(node_id) == low_deg
    //         [Tij1j2, Tij1j2_tilde] = ...
    //             make_Tij1j2s_edges( ...
    //             node_id, T_edges, Tijs_scaled, edges, problem_data);
    //         Tij_2deg_recovery = cat(3, Tij_2deg_recovery, Tij1j2);
    //         Tij_tilde_2deg_recovery = cat( ...
    //             3, Tij_tilde_2deg_recovery, Tij1j2_tilde);
    int lowDegIdx = 0;
    for (int nodeId = 0; nodeId < n; ++nodeId)
    {
        if (nodeDegrees(nodeId, 0) == lowDeg)
        {
            SomUtils::MatD Tij1j2(SomUtils::MatD::Zero(d, lowDeg));
            SomUtils::MatD Tij1j2Tilde(SomUtils::MatD::Zero(p, lowDeg));
            Prob.makeTij1j2sEdges(nodeId, nodeDegrees, Tedges, Tij1j2, Tij1j2Tilde);
            Tij2degRecovery[lowDegIdx] = Tij1j2;
            TijTilde2degRecovery[lowDegIdx] = Tij1j2Tilde;
            lowDegIdx++;
        }
    }

    // Tij_tilde_2deg_recovery=multiprod(Qalign, Tij_tilde_2deg_recovery);
    auto TijTilde2degRecoveryQalign = TijTilde2degRecovery;
    SomUtils::VecMatD RmanoptOutQalign(lowDeg, SomUtils::MatD::Zero(p, d));
    lowDegIdx = 0;
    for (int nodeId = 0; nodeId < n; ++nodeId)
    {
        if (nodeDegrees(nodeId, 0) == lowDeg)
        {
            TijTilde2degRecoveryQalign[lowDegIdx] = Qalign * TijTilde2degRecovery[lowDegIdx];
            RmanoptOutQalign[lowDegIdx] = Qalign * RmanoptOut[nodeId];
            lowDegIdx++;
        }
    }
    // RitildeEst = RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
    SomUtils::VecMatD RitildeEst(numNodesLowDeg, SomUtils::MatD::Zero(p, d));
    SomUtils::VecMatD Qxs(numNodesLowDeg, SomUtils::MatD::Zero(p, p));
    SomUtils::VecMatD Qbs(numNodesLowDeg, SomUtils::MatD::Zero(p, p));
    // ROFL_VAR3(RmanoptOutQalign.size(), TijTilde2degRecoveryQalign.size(), RitildeEst.size())
    // ROFL_VAR2(Qxs.size(), Qbs.size())
    // ROFL_VAR2(RmanoptOutQalign[0].rows(), RmanoptOutQalign[0].cols())
    // ROFL_VAR2(RmanoptOutQalign[1].rows(), RmanoptOutQalign[1].cols())
    // ROFL_VAR2(TijTilde2degRecoveryQalign[0].rows(), TijTilde2degRecoveryQalign[0].cols())
    // ROFL_VAR2(TijTilde2degRecoveryQalign[1].rows(), TijTilde2degRecoveryQalign[1].cols())
    // ROFL_VAR2(RitildeEst[0].rows(), RitildeEst[0].cols())
    // ROFL_VAR2(RitildeEst[1].rows(), RitildeEst[1].cols())

    // return 0;
    Prob.RbRecovery(RmanoptOutQalign, TijTilde2degRecoveryQalign, RitildeEst, Qxs, Qbs);
    ROFL_VAR1("RitildeEst")
    for (int i = 0; i < numNodesLowDeg; ++i)
    {
        ROFL_VAR1(RitildeEst[i])
    }

    // R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);
    SomUtils::VecMatD Rrecovered(n, SomUtils::MatD::Zero(d, d));
    lowDegIdx = 0;
    for (int nodeId = 0; nodeId < n; ++nodeId)
    {
        if (nodeDegrees(nodeId, 0) == lowDeg)
        {
            Rrecovered[nodeId] = RitildeEst[lowDegIdx].block(0, 0, d, d);
            lowDegIdx++;
        }
        else
        {
            // R_tilde2_edges = multiprod(repmat(Qalign, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));
            auto tmp = Qalign * RmanoptOut[nodeId];
            Rrecovered[nodeId] = tmp.block(0, 0, d, d);
        }
    }

    // low_deg_nodes_ids = find(problem_data.node_degrees <= low_deg); %[1 5]'
    // for ii = 1:N
    //     if ismember(ii, low_deg_nodes_ids)
    lowDegIdx = 0;
    for (int nodeId = 0; nodeId < n; ++nodeId)
    {
        ROFL_VAR2(nodeId, Rrecovered[nodeId].determinant())

        if (nodeDegrees(nodeId, 0) == lowDeg)
        {
            // id_low_deg = find(low_deg_nodes_ids == ii);
            // P_i = recover_R_deg2(Tij_tilde_2deg_recovery, id_low_deg, d);
            // R_recovered(:,:,ii) = P_i * R_recovered(:,:,ii);

            if (Rrecovered[nodeId].determinant() < 0)
            {
                SomUtils::MatD Pi = SomUtils::MatD::Zero(3, 3);
                Prob.recoverRdeg2(TijTilde2degRecovery, lowDegIdx, Pi);
                Rrecovered[nodeId] = Pi * Rrecovered[nodeId];
            }
            lowDegIdx++;
        }
    }

    // disp("multidet(R_recovered)")
    // disp(multidet(R_recovered))
    for (int i = 0; i < n; ++i)
    {
        ROFL_VAR3(i, Rrecovered[i], Rrecovered[i].determinant())
    }

    // R_recovered(:,:,nodes_high_deg) = R_tilde2_edges(1:d,:,:);
    // nodes_low_deg = ~nodes_high_deg;

    // T_diffs_shifted = Qalign * T_edges; %this has last row to 0
    SomUtils::MatD TedgesShifted(SomUtils::MatD::Zero(p, numEdges));
    for (int i = 0; i < numEdges; ++i)
    {
        TedgesShifted.col(i) = Qalign * Tedges.col(i);
    }
    ROFL_VAR2(TedgesShifted, TedgesShifted.block(3, 0, p - 3, numEdges).cwiseAbs().maxCoeff());

    SomUtils::MatD Trecovered(SomUtils::MatD::Zero(p, n));
    Prob.edgeDiffs2T(0, TedgesShifted, n, Trecovered);

    SomUtils::MatD LambdasRecovered(SomUtils::MatD::Zero(numEdges, 1));
    double lambdasCoeff = LambdasGtEig(0, 0) / LambdasManoptOut(0, 0);
    LambdasRecovered = lambdasCoeff * LambdasManoptOut;

    ROFL_VAR2(LambdasRecovered.transpose(), LambdasGtEig.transpose())

    ROFL_VAR1(Prob.costEigen(RgtEig, TgtEig, LambdasGtEig))

    ROFL_VAR1(Prob.costEigen(RmanoptOut, TmanoptOut, LambdasManoptOut))

    ROFL_VAR1(Prob.costEigen(Rrecovered, Trecovered.block(0, 0, d, n), LambdasRecovered))

    // T_recovered_pre = recover_T_edges(T_diffs_shifted(1:d,:), ...
    //     edges, d, problem_data.node_degrees, low_deg, Tij_tilde_2deg_recovery);
    // T_recovered = edge_diffs_2_T(T_recovered_pre, edges, N);
    // % T_recovered = edge_diffs_2_T(T_diffs_shifted(1:d, :), edges, N);
    // %%

    // X_recovered.R = R_recovered;
    // X_recovered.T = T_recovered;
    // X_recovered.lambda = lambdas_recovered;

    // problem_data_next = problem_data; %TODO: fix this line after recovery works
    // cost_out = ssom_cost(X_recovered, problem_data_next);
    // disp("cost_out AFTER RECOVERY")
    // disp(cost_out)

    return 0;
}

procrustesOutput procrustesCodemeta(const SomUtils::MatD &X, const SomUtils::MatD &Y, double &d, SomUtils::MatD &Z, bool doScaling, std::string &doReflection)
{
    // bool doScaling=true;
    // doReflection='best';
    // %optional parameters
    // ivarargin=1;
    // while(ivarargin<=length(varargin))
    //     switch(lower(varargin{ivarargin}))
    //         case 'scaling'
    //             ivarargin=ivarargin+1;
    //             doScaling=varargin{ivarargin};
    //         case 'reflection'
    //             ivarargin=ivarargin+1;
    //             doReflection=varargin{ivarargin};
    //         otherwise
    //             error(['Argument ' varargin{ivarargin} ' not valid!'])
    //     end
    //     ivarargin=ivarargin+1;
    // end

    // [n, m]   = size(X);
    // [ny, my] = size(Y);
    int xrows = X.rows();
    int xcols = X.cols();
    int yrows = Y.rows();
    int ycols = Y.cols();

    // Check input sizes.
    if (xrows <= 0 || xcols <= 0 || yrows <= 0 || ycols <= 0)
        throw std::invalid_argument("procrustes:InputSizeMismatch: X and Y must be non-empty matrices.");

    // if ny ~= n
    //     error('procrustes:InputSizeMismatch',...
    //         'X and Y must have the same number of rows (points).');
    // elseif my > m
    //     error('procrustes:InputSizeMismatch',...
    //         'Y cannot have more columns (variables) than X.');
    // end
    if (yrows != xrows)
        throw std::invalid_argument("procrustes:InputSizeMismatch: X and Y must have the same number of rows (points).");

    if (ycols > xcols)
        throw std::invalid_argument("procrustes:InputSizeMismatch: Y cannot have more columns (variables) than X.");

    // % Center at the origin.
    // muX = mean(X,1);
    // muY = mean(Y,1);
    // X0 = X - repmat(muX, n, 1);
    // Y0 = Y - repmat(muY, n, 1);
    SomUtils::MatD muX = X.colwise().mean();
    SomUtils::MatD muY = Y.colwise().mean();
    SomUtils::MatD X0 = X;
    SomUtils::MatD Y0 = Y;
    for (int i = 0; i < xrows; ++i)
    {
        X0.row(i) -= muX;
        Y0.row(i) -= muY;
    }
    // ssqX = sum(X0.^2,1);
    // ssqY = sum(Y0.^2,1);
    // constX = all(ssqX <= abs(eps(class(X))*n*muX).^2);
    // constY = all(ssqY <= abs(eps(class(X))*n*muY).^2);
    // ssqX = sum(ssqX);
    // ssqY = sum(ssqY);
    double ssqX = X0.squaredNorm();
    double ssqY = Y0.squaredNorm();
    bool constX = (ssqX <= std::numeric_limits<double>::epsilon() * xrows * muX.array().abs().square().sum());
    bool constY = (ssqY <= std::numeric_limits<double>::epsilon() * xrows * muY.array().abs().square().sum());
    // HP) ssqX = sum(ssqX), ssqY = sum(ssqY) are never needed
    //////////////////////////////////////////////////////////////

    // if ~constX & ~constY
    if (!constX && !constY)
    {
        //     % The "centered" Frobenius norm.
        //     normX = sqrt(ssqX); % == sqrt(trace(X0*X0'))
        //     normY = sqrt(ssqY); % == sqrt(trace(Y0*Y0'))
        double normX = std::sqrt(ssqX); // == sqrt(trace(X0*X0'))
        double normY = std::sqrt(ssqY); // == sqrt(trace(Y0*Y0'))

        // % Scale to equal (unit) norm.
        // X0 = X0 / normX;
        // Y0 = Y0 / normY;
        X0 /= normX;
        Y0 /= normY;

        // % If Y has fewer variables than X, pad it with zeros.
        // % Make sure they're in the same dimension space.
        // if my < m
        //     Y0 = [Y0 zeros(n, m-my)];
        if (ycols < xcols)
            Y0.conservativeResize(Eigen::NoChange, xcols);

        // % The optimum rotation matrix of Y.
        // A = X0' * Y0;
        SomUtils::MatD A = X0.transpose() * Y0;
        // [L, D, M] = svd(A);
        Eigen::JacobiSVD<SomUtils::MatD> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
        SomUtils::MatD L = svd.matrixU();                     // Left singular vectors
        SomUtils::MatD D = svd.singularValues().asDiagonal(); // Singular values
        SomUtils::MatD M = svd.matrixV();                     // Right singular vectors
        // T = M * L';
        SomUtils::MatD T = M * L.transpose(); // Rotation matrix

        // if isempty(doReflection) % 'best'
        //     % Let the data decide if a reflection is needed.
        // else
        //     haveReflection = (det(T) < 0);
        //     % If we don't have what was asked for ...
        //     if (doReflection ~= haveReflection)
        //         % ... then either force a reflection, or undo one.
        //         M(:,end) = -M(:,end);
        //         D(end,end) = -D(end,end);
        //         T = M * L';
        if (doReflection.empty())
        {
            // Let the data decide if a reflection is needed.
        }
        else
        {
            bool haveReflection = (T.determinant() < 0);
            if (doReflection != "haveReflection")
            {
                M.col(xcols - 1) = -M.col(xcols - 1);
                D(xcols - 1, xcols - 1) = -D(xcols - 1, xcols - 1);
                T = M * L.transpose();
            }
        }

        // % The minimized unstandardized distance D(X0,b*Y0*T) is
        // % ||X0||^2 + b^2*||Y0||^2 - 2*b*trace(T*X0'*Y0)
        // traceTA = sum(diag(D)); % == trace(sqrtm(A'*A)) when doReflection is 'best'
        double traceTA = D.diagonal().sum();

        // if doScaling
        if (doScaling)
        {
            //     % The optimum scaling of Y.
            //     b = traceTA * normX / normY;
            double b = traceTA * normX / normY;
            //     % The standardized distance between X and b*Y*T+c.
            //     d = 1 - traceTA.^2;
            double d = 1 - traceTA * traceTA;
            //     if nargout > 1
            //         Z = normX*traceTA * Y0 * T + repmat(muX, n, 1);
            SomUtils::MatD Z = normX * traceTA * Y0 * T;
            for (int i = 0; i < Z.rows(); ++i)
            {
                Z.row(i) += muX;
            }
            //     else % if ~doScaling
            //         b = 1;
            //         % The standardized distance between X and Y*T+c.
            //         d = 1 + ssqY/ssqX - 2*traceTA*normY/normX;
            //         if nargout > 1
            //             Z = normY*Y0 * T + repmat(muX, n, 1);
            if (!doScaling)
            {
                b = 1;
                // The standardized distance between X and Y*T+c.
                double d = 1 + ssqY / ssqX - 2 * traceTA * normY / normX;
                Z = normY * Y0 * T;
                for (int i = 0; i < Z.rows(); ++i)
                {
                    Z.row(i) += muX;
                }
            }

            // if nargout > 2
            //     if my < m
            //         T = T(1:my,:);
            //     c = muX - b*muY*T;
            //     transform = struct('T',T, 'b',b, 'c',repmat(c, n, 1));
            if (ycols < xcols)
            {
                auto tmp = T.block(0, 0, ycols, xcols);
                T.resize(ycols, xcols);
                T = tmp;
            }
            SomUtils::MatD c = muX - b * muY * T;
            SomUtils::MatD Zc = SomUtils::MatD::Zero(xrows, xcols);
            for (int i = 0; i < xrows; ++i)
            {
                Zc.row(i) = c.transpose();
            }
            T = Zc;

            // transform = struct('T',T, 'b',b, 'c',Zc);
        }
        return procrustesOutput(T, d, Z);
    }
    else if (constX) //////////////////////////////////////////////////////////////
    {
        // % The degenerate cases: X all the same, and Y all the same.
        // elseif constX
        //     d = 0;
        //     Z = repmat(muX, n, 1);
        //     T = eye(my,m);
        //     transform = struct('T',T, 'b',0, 'c',Z);
        double d = 0;
        SomUtils::MatD T = SomUtils::MatD::Identity(ycols, xcols);
        SomUtils::MatD Z = SomUtils::MatD::Zero(xrows, xcols);
        for (int i = 0; i < xrows; ++i)
        {
            Z.row(i) = muX;
        }
        return procrustesOutput(T, d, Z);
    }
    else
    {
        // else % ~constX & constY
        //     d = 1;
        //     Z = repmat(muX, n, 1);
        //     T = eye(my,m);
        //     transform = struct('T',T, 'b',0, 'c',Z);
        double d = 1;
        SomUtils::MatD T = SomUtils::MatD::Identity(ycols, xcols);
        SomUtils::MatD Z = SomUtils::MatD::Zero(xrows, xcols);
        for (int i = 0; i < xrows; ++i)
        {
            Z.row(i) = muX;
        }
        return procrustesOutput(T, d, Z);
    }

    return procrustesOutput(); // should never be reached
}

void devectorizeXeig(const SomUtils::MatD &xEig, SomUtils::VecMatD &R, SomUtils::MatD &T, SomUtils::MatD &Lambdas, int p, int d, int n, int numEdges)
{
    ROFL_ASSERT(xEig.size() == p * d * n + p * n + numEdges)

    SomUtils::MatD RhSt = SomUtils::MatD::Zero(p, d * n);
    RhSt = xEig.block(0, 0, p * d * n, 1).reshaped(p, d * n);
    SomUtils::unStackH(RhSt, R, d);
    T = xEig.block(p * d * n, 0, p * n, 1).reshaped(p, n);
    Lambdas = xEig.block(p * d * n + p * n, 0, numEdges, 1);
}