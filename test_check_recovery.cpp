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

    // load('poc2degree_data/recovery_check.mat','R');
    // load('poc2degree_data/recovery_check.mat','T');
    // R_manopt_out = R;
    // T_manopt_out = T;
    // X_manopt_out.R = R_manopt_out;
    // X_manopt_out.T = T_manopt_out;
    // clear R T;
    // load('poc2degree_data/recovery_check.mat','params');
    // load('poc2degree_data/recovery_check.mat','edges');
    // load('poc2degree_data/recovery_check.mat','d');
    // load('poc2degree_data/recovery_check.mat','N');
    // load('poc2degree_data/recovery_check.mat','staircase_step_idx');
    // load('poc2degree_data/recovery_check.mat','Tijs');
    // load('poc2degree_data/recovery_check.mat','X_gt');
    // load('poc2degree_data/recovery_check.mat','problem_struct_next');

    int staircaseStepIdx = 5;
    int nrs = staircaseStepIdx - 1;
    int d = 3;
    int n = 6;

    int numEdges = 18;

    SomUtils::SomSize somSzD(d, d, n);
    SomUtils::SomSize somSzNrs(nrs, d, n);
    SomUtils::MatD Tijs(d, numEdges);
    Eigen::MatrixXi edges(numEdges, 2);

    SomUtils::readCsvTijs("../data/recov_tijs.csv", Tijs, d, numEdges);
    SomUtils::readCsvEdges("../data/recov_edges.csv", edges);

    // problem d x d x n
    integer numoftypes = 2; // 2 i.e. (3D) Stiefel + Euclidean
    integer numofmani1 = n; // num of Stiefel manifolds
    integer numofmani2 = 1;
    ROPTLIB::Stiefel mani1(d, d);
    mani1.ChooseParamsSet2();
    ROPTLIB::Euclidean mani2(d, n);
    ROPTLIB::ProductManifold ProdMani(numoftypes, &mani1, numofmani1, &mani2, numofmani2);
    // Read from csv
    ROPTLIB::Vector xGt = ProdMani.RandominManifold();
    SomUtils::readCsvInitguess("../data/recov_x_gt.csv", xGt);
    xGt.Print("xGt");
    ROPTLIB::SampleSomProblem Prob(somSzD, Tijs, edges);
    Prob.SetDomain(&ProdMani);

    // ROPT to Eig (GT)
    SomUtils::MatD XgtVecEig(SomUtils::MatD::Zero(d * d * n + d * n, 1));
    SomUtils::VecMatD RgtEig(n, SomUtils::MatD::Zero(d, d));
    SomUtils::MatD TgtEig(SomUtils::MatD::Zero(d, n));
    Prob.RoptToEig(xGt, XgtVecEig);
    Prob.getRotations(XgtVecEig, RgtEig);
    Prob.getTranslations(XgtVecEig, TgtEig);
    Prob.setGt(RgtEig, TgtEig);

    for (auto &m : RgtEig)
        ROFL_VAR1(m)
    ROFL_VAR1(TgtEig)

    // problem nrs x d x n
    ROPTLIB::Stiefel mani1nrs(nrs, d);
    mani1.ChooseParamsSet2();
    ROPTLIB::Euclidean mani2nrs(nrs, n);
    ROPTLIB::ProductManifold ProdManiNrs(numoftypes, &mani1nrs, numofmani1, &mani2nrs, numofmani2);
    // Read from csv
    ROPTLIB::Vector xManoptOut = ProdManiNrs.RandominManifold();
    SomUtils::readCsvInitguess("../data/recov_x_manopt_out.csv", xManoptOut);
    xManoptOut.Print("xManoptOut");
    ROPTLIB::SampleSomProblem ProbNrs(somSzNrs, Tijs, edges);
    ProbNrs.SetDomain(&ProdManiNrs);

    // ROPT to Eig (Manopt Output)
    SomUtils::MatD XmanoptOutVecEig(SomUtils::MatD::Zero(nrs * d * n + nrs * n, 1));
    SomUtils::VecMatD RmanoptOutEig(n, SomUtils::MatD::Zero(nrs, d));
    SomUtils::MatD TmanoptOutEig(SomUtils::MatD::Zero(nrs, n));
    ProbNrs.RoptToEig(xManoptOut, XmanoptOutVecEig);
    ProbNrs.getRotations(XmanoptOutVecEig, RmanoptOutEig);
    ProbNrs.getTranslations(XmanoptOutVecEig, TmanoptOutEig);
    Prob.setGt(RgtEig, TgtEig); // !!

    for (auto &m : RmanoptOutEig)
        ROFL_VAR1(m)
    ROFL_VAR1(TmanoptOutEig)

    // back to SE(d)^N
    SomUtils::VecMatD Rrecovered(n, SomUtils::MatD::Zero(d, d));
    SomUtils::MatD Trecovered(SomUtils::MatD::Zero(d, n));
    bool recSEDNsuccess = ProbNrs.recoverySEdN(staircaseStepIdx,
                                               RmanoptOutEig, TmanoptOutEig,
                                               Rrecovered, Trecovered);

    ROFL_VAR1("Printing R, T SE(d)^N")
    for (auto &m : Rrecovered)
        ROFL_VAR1(m)
    ROFL_VAR1(Trecovered)

    ROFL_VAR1(recSEDNsuccess)

    // globalize
    SomUtils::VecMatD Rout(n, SomUtils::MatD::Zero(d, d));
    SomUtils::MatD Tout(SomUtils::MatD::Zero(d, n));
    int src = 0;
    bool globalRecoverySuccess = ProbNrs.globalize(src, Rrecovered, Trecovered,
                                                   Rout, Tout);

    ROFL_VAR1("Printing R, T out")
    for (auto &m : Rout)
        ROFL_VAR1(m)
    ROFL_VAR1(Tout)

    ROFL_VAR1(globalRecoverySuccess)

    // % if staircase_step_idx > d+1

    // low_deg = 2; %TODO: not necessarily in more complex graph cases
    // nodes_high_deg = params.node_degrees > low_deg;

    // [T_edges, T1_offset] = make_T_edges(T_manopt_out, edges);

    // RT_stacked_high_deg = [matStackH(R_manopt_out(:,:,nodes_high_deg)), T_edges];
    // Qx_edges = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);
    // % RT_stacked_high_deg_poc = Qx_edges * RT_stacked_high_deg;

    // R_tilde2_edges = multiprod(repmat(Qx_edges, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));

    // R_recovered = zeros(d,d,N);
    // R_recovered(:,:,nodes_high_deg) = R_tilde2_edges(1:d,:,:);

    // %% check_Tij1j2_init on node 1

    // node_id = 1;
    // Ritilde2 = R_manopt_out(:,:,node_id);

    // check_Rb_params.node_degrees = params.node_degrees;
    // check_Rb_params.nrs = staircase_step_idx - 1;
    // check_Rb_params.d = d;
    // check_Rb_params.R_gt = params.R_gt;
    // check_Rb_params.T_gt = params.T_gt;

    // [Tij1j2, Tij1j2_tilde] = make_Tij1j2s_edges(node_id, Qx_edges * T_edges, Tijs, edges,check_Rb_params);

    // [~,~,~,RbEst]=recoverRitilde(Ritilde2,Tij1j2_tilde);
    // Qb = blkdiag(eye(2),RbEst');

    // Qx=align2d(Tij1j2_tilde);
    // check_Tij1j2_init(1, multiprod(Qx_edges, R_manopt_out), Tij1j2, Tij1j2_tilde, ...
    //     Qx, Qb);

    // %% run only later

    // % nodes_low_deg = ~nodes_high_deg;
    // for node_id = 1:length(params.node_degrees)
    //     node_deg = params.node_degrees(node_id);
    //     if node_deg == low_deg
    //         fprintf("Running recoverRitilde() on node %g\n", node_id);
    //         R_i_tilde2 = R_manopt_out(:,:,node_id);

    //         check_Rb_params.node_degrees = params.node_degrees;
    //         check_Rb_params.nrs = staircase_step_idx - 1;
    //         check_Rb_params.d = d;
    //         check_Rb_params.R_gt = params.R_gt;
    //         check_Rb_params.T_gt = params.T_gt;

    //         cost_gt = rsom_cost_base(X_gt, problem_struct_next);
    //         disp("cost_gt")
    //         disp(cost_gt)

    //         cost_manopt_output = rsom_cost_base(X_manopt_out, problem_struct_next);
    //         disp("cost_manopt_output")
    //         disp(cost_manopt_output)

    //         T_diffs_shifted = Qx_edges * T_edges; %this has last row to 0
    //         [Tij1j2, Tij1j2_tilde] = make_Tij1j2s_edges(node_id, T_diffs_shifted, Tijs,edges,params);

    //         [RitildeEst1,RitildeEst2,Qx,Rb_i] = ...
    //             recoverRitilde(Qx_edges* R_i_tilde2,Tij1j2_tilde);
    //         disp('')
    //         % TODO: how to decide between RitildeEst1,RitildeEst2??
    //         if (det(RitildeEst1(1:d,:)) > 0)
    //             R_recovered(:,:,node_id) = RitildeEst1(1:d,:);
    //         else
    //             R_recovered(:,:,node_id) = RitildeEst2(1:d,:);
    //         end

    //         disp("Checking make_T_edges <-> edge_diffs_2_T");

    //         T_recovered = edge_diffs_2_T(T_diffs_shifted(1:d,:), edges, N);
    //         disp('')
    //     end
    // end

    // %checking that cost has not changed during "recovery"
    // X_recovered.T = T_recovered;
    // X_recovered.R = R_recovered;
    // cost_out = rsom_cost_base(X_recovered, problem_struct_next);
    // disp("cost_out")
    // disp(cost_out)

    // % transf_out = RT2G(R_recovered, T_recovered); %rsom_genproc() function output

    // disp("[matStackH(X_gt.R); matStackH(R_recovered)]");
    // disp([matStackH(X_gt.R); matStackH(R_recovered)]);

    // R_global = R_recovered(:,:,1) * X_gt.R(:,:,1)'; %!!
    // % code for making all rotations global at once
    // R_recovered_global = multiprod(repmat(R_global', 1, 1, N), R_recovered);
    // disp("[matStackH(X_gt.R); matStackH(R_recovered_global)]");
    // disp([matStackH(X_gt.R); matStackH(R_recovered_global)]);

    // T_global = R_global * T_recovered(:,1) - X_gt.T(:,1); %!!
    // % code for making all translation global at once
    // disp("[X_gt.T; T_recovered]");
    // T_recovered_global = R_global' * T_recovered - T_global;
    // disp([X_gt.T; T_recovered_global]);

    // for ii = 1:N
    //     R_gt_i = X_gt.R(:,:,ii);
    //     R_recov_i_global = R_recovered_global(:,:,ii); %GLOBAL!
    //     fprintf("ii %g\n", ii);
    //     % rotations
    //     disp("R_gt_i, R_recov_i");
    //     disp([R_gt_i, R_recov_i_global]);
    //     disp("is_equal_floats(R_gt_i, R_recov_i_global)")
    //     disp(is_equal_floats(R_gt_i, R_recov_i_global))
    //     if (~is_equal_floats(R_gt_i, R_recov_i_global))
    //         error("rot found NOT equal")
    //     end
    //     % translations
    //     T_gt_i = X_gt.T(:,ii);
    //     T_recov_i_global = T_recovered_global(:,ii);
    //     disp("[X_gt.T, T_recovered]");
    //     disp([T_gt_i, T_recov_i_global]);
    //     disp("is_equal_floats(T_gt_i, T_recov_i_global)")
    //     disp(is_equal_floats(T_gt_i, T_recov_i_global))
    //     if (~is_equal_floats(T_gt_i, T_recov_i_global))
    //         error("transl found NOT equal")

    //     end
    // end

    // end %file function
}