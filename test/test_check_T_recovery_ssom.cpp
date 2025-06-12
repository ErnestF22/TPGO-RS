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

void devectorizeXeig(const SomUtils::MatD &xEig, SomUtils::VecMatD &R, SomUtils::MatD &T, SomUtils::MatD &Lambdas, int p, int d, int n, int numEdges)
{
    ROFL_ASSERT(xEig.size() == p * d * n + p * n + numEdges)

    SomUtils::MatD RhSt = SomUtils::MatD::Zero(p, d * n);
    RhSt = xEig.block(0, 0, p * d * n, 1).reshaped(p, d * n);
    SomUtils::unStackH(RhSt, R, d);
    T = xEig.block(p * d * n, 0, p * n, 1).reshaped(p, n);
    Lambdas = xEig.block(p * d * n + p * n, 0, numEdges, 1);
}

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
    SomUtils::MatD Qalign(SomUtils::MatD::Zero(p, p));
    Prob.align3d(RTstackedHighDeg, Qalign);
    ROFL_VAR1(Qalign)

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
    auto RmanoptOutQalign = RmanoptOut;
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
    Prob.RbRecovery(RmanoptOutQalign, TijTilde2degRecoveryQalign, RitildeEst, Qxs, Qbs);
    ROFL_VAR1("RitildeEst")
    for (int i = 0; i < numNodesLowDeg; ++i)
    {
        ROFL_VAR1(RitildeEst[i].transpose())
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
            auto tmp = Qalign * RmanoptOut[nodeId];
            Rrecovered[nodeId] = tmp.block(0, 0, d, d);
        }
    }

    // R_tilde2_edges = multiprod(repmat(Qalign, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));

    // low_deg_nodes_ids = find(problem_data.node_degrees <= low_deg); %[1 5]'
    // for ii = 1:N
    //     if ismember(ii, low_deg_nodes_ids)
    //         id_low_deg = find(low_deg_nodes_ids == ii);
    //         P_i = recover_R_deg2(Tij_tilde_2deg_recovery, id_low_deg, d);
    //         R_recovered(:,:,ii) = P_i * R_recovered(:,:,ii);
    //     end
    // end

    // R_recovered(:,:,nodes_high_deg) = R_tilde2_edges(1:d,:,:);
    // nodes_low_deg = ~nodes_high_deg;

    // disp("multidet(R_recovered)")
    // disp(multidet(R_recovered))

    // % T_edges_recovered = recover_T_edges(Qalign * T_edges, edges, ...
    // %     node_degrees, low_deg, Qxs, Qbs, Qalign, Qx_edges);
    // T_diffs_shifted = Qalign * T_edges; %this has last row to 0
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
