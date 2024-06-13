function test_full_recovery


N = 6;
load('Qbnn_data/testdata.mat','edges')
adj_mat = make_adj_mat_from_edges(edges, N);
testdata = testNetwork_adj(3, adj_mat, 'banded', 3);
G = graph(make_adj_mat_from_edges(testdata.E,N));
figure(1000)
plot(G)
title('graph')


load('Qbnn_data/testdata.mat','R')
load('Qbnn_data/testdata.mat','T')
load('Qbnn_data/testdata.mat','Tijs')

nrs = size(R, 1);



d = size(R, 2);

node_degrees = [2,2,4,3,3,4];
nodes_high_deg = node_degrees == 3 | node_degrees == 4;
nodes_low_deg = ~nodes_high_deg;

% R_stacked_high_deg = matStackH(R(:,:,nodes_high_deg));
% Qx = POCRotateToMinimizeLastEntries(R_stacked_high_deg);
% R_stacked_high_deg_poc = Qx * R_stacked_high_deg;
% R_high_deg_poc = matUnstackH(R_stacked_high_deg_poc, d);
%
% disp('R_stacked_high_deg_poc')
% disp(R_stacked_high_deg_poc)

% disp("max(R_stacked_high_deg_poc(d+1:end, :), [], ""all"")")
% disp(max(R_stacked_high_deg_poc(d+1:end, :), [], "all"))



T_edges = make_T_edges(T, edges);

RT_stacked_high_deg = [matStackH(R(:,:,nodes_high_deg)), T_edges];
Qx = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);
RT_stacked_high_deg_poc = Qx * RT_stacked_high_deg;

% R_recovered = zeros(nrs,d,N);
% R_recovered(:,:,nodes_high_deg) = R_high_deg_poc;
% R_recovered(:,:,nodes_low_deg) = multiprod(repmat(Qx, 1, 1, sum(nodes_low_deg)), R(:,:,nodes_low_deg));
R_recovered = multiprod(repmat(Qx, 1, 1, N), R);


disp('RT_stacked_high_deg_poc')
disp(RT_stacked_high_deg_poc)

disp("max(RT_stacked_high_deg_poc(d+1:end, :), [], ""all"")")
disp(max(RT_stacked_high_deg_poc(d+1:end, :), [], "all"))



X.R = R;
X.T = T;
problem_data.nrs = 4;
problem_data.Tijs = Tijs;
problem_data.d = 3;
problem_data.N = 6;
problem_data.edges = edges;
disp("rsom_cost_base(X, problem_data)")
disp(rsom_cost_base(X, problem_data))


for e = 1:size(T_edges, 2)
    ii = edges(e,1);
    %     if nodes_low_deg(ii)
    %         continue;
    %     end
    jj = edges(e,2);
    %
    disp([ii, jj])
    Tij = Tijs(:,e);
    %     disp(R(:,:,ii) * Tijs(:, e) - T(:,jj) + T(:,ii))
    R_i = R_recovered(:,:,ii);
    %     R_j = R_high_deg_poc(:,:,jj);
    RiTij = R_i * Tij;
    disp(RiTij)
    if(RiTij(end) > 1e-5)
        T_e = T_edges(:, e);
        disp(R_i * Tij + T_e)
        disp(Qx * R_i * Tij + Qx * T_e)
        disp("error")
    end
    %
end


sz = [nrs, d, N];
[R_final_stiefel, ~] = ...
    Qc_recovery_Rb_initguess(sz, edges, R, T, Tijs, node_degrees);
R_final_stiefel(:,:,nodes_high_deg) = R_recovered(:,:,nodes_high_deg);

load('/home/ernest/workspace/matlab_ws/shape_of_motion/genproc/Qbnn_data/R_gt.mat', 'R_gt')
R_final = R_final_stiefel(1:d, :, :);

for ii = 1:N
    disp('ii')
    disp(ii)
    disp("R_final(:,:,ii) * R_gt(:,:,ii)'")
    disp(R_final(:,:,ii) * R_gt(:,:,ii)')
end

end %file function
