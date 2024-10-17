low_deg = 2;

%computing Tij12

N = 6;
load('poc2degree_data/ws.mat','edges')
adj_mat = make_adj_mat_from_edges(edges, N);
testdata = testNetwork_adj(3, adj_mat, 'banded', 3);
G = graph(make_adj_mat_from_edges(testdata.E,N));
figure(1000)
plot(G)
title('graph')


load('poc2degree_data/ws.mat','R')
load('poc2degree_data/ws.mat','T')
load('poc2degree_data/ws.mat','Tijs')

nrs = size(R, 1);
d = size(R, 2);

node_degrees = [2,2,4,3,3,4];
nodes_high_deg = node_degrees == 3 | node_degrees == 4;
nodes_low_deg = ~nodes_high_deg;



%%

X.R = R;
X.T = T;
problem_data.nrs = 4;
problem_data.Tijs = Tijs;
problem_data.d = 3;
problem_data.N = 6;
problem_data.edges = edges;
disp("rsom_cost_base(X, problem_data)")
disp(rsom_cost_base(X, problem_data))


load('poc2degree_data/R_gt.mat', "R_globalframe")
load('poc2degree_data/T_gt.mat', "T_globalframe")
R_gt = R_globalframe;
T_gt = T_globalframe;
X_gt.R = R_gt;
X_gt.T = T_gt;
% load('poc2degree_data/problem_data_gt.mat', "problem_data_gt")
disp("rsom_cost_base(X_gt, problem_data)")
disp(rsom_cost_base(X_gt, problem_data))


%% grid search test


angles_rad = 0:1:360;
rot_dist_best = 1e+6;
R_gt_stiefel = [R_gt(:,:,2); zeros(1,3)];

Qxs = randrot_manopt(4, 1000);

for a = angles_rad
    Rb = rot2d(pi*a/180.0);
    Qb = blkdiag(eye(2), Rb);
    for jj = 1:size(Qxs,3)
        Qx = Qxs(:,:,jj);
        Ri2_tilde_cand = Qx' * Qb * Qx * R_gt_stiefel;
        rot_dist_i = norm(R(:,:,2)-Ri2_tilde_cand,'fro');
        if (rot_dist_i < rot_dist_best)
            rot_dist_best = rot_dist_i;
            Qb_best = Qb;
            Qx_best = Qx;
            disp("Found better approximation")
            disp(norm(R(:,:,2)-Ri2_tilde_cand,'fro'));
        end
    end
end

pinv(R_gt_stiefel) * R(:,:,2)


function m = rot2d(angle_rad)
    m = [cos(angle_rad), -sin(angle_rad); sin(angle_rad), cos(angle_rad)];
end

function m_out = so4_sampling(A)
%         A = randn(n);
        [Q, RR] = qr(A);
        Q = Q * diag(sign(diag(RR))); %% Mezzadri 2007
        
        % If Q is in O(n) but not in SO(n), we permute the two first
        % columns of Q such that det(new Q) = -det(Q), hence the new Q will
        % be in SO(n), uniformly distributed.
        if det(Q) < 0
            Q(:, [1 2]) = Q(:, [2 1]);
        end
        
        m_out = Q;
end