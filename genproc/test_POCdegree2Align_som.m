function test_POCdegree2Align_som
% From work with Ernesto
% See equations (41)-(45) in our working document
% flagRandomizeQx=true;
% flagPerturbMeasurements=false;

resetRands()


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




%%
T_edges = make_T_edges(T, edges);
RT_stacked_high_deg = [matStackH(R(:,:,nodes_high_deg)), T_edges];
Qx = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);
RT_stacked_high_deg_poc = Qx * RT_stacked_high_deg;
disp("RT_stacked_high_deg_poc")
disp(RT_stacked_high_deg_poc)
disp("nodes low deg")
disp(nodes_low_deg)
col_ids = reshape(1:12, d, []);
for ii = 1:4
    disp(ii)
    %     A = R_gt(:,:,ii+2) * RT_stacked_high_deg_poc(1:d, col_ids(:,ii))'; %works only in this particular case
    A = RT_stacked_high_deg_poc(1:d, col_ids(:,ii)) * R_gt(:,:,ii+2)'; %works only in this particular case
    disp(A);
end

%%

Qs_poc = [];

for deg_i = 1:length(node_degrees)
    ii = node_degrees(deg_i);
    if ii == low_deg
        fprintf("Running recoverRitilde() on node %g\n", deg_i);
        R_i = R(:,:,deg_i);

        %computing Tij12
        Tij1j2 = zeros(d,low_deg);
        num_js_found = 0;
        for e = 1:size(edges)
            e_i = edges(e,1);
            %         e_j = edges(e,2);
            if (e_i == ii)
                num_js_found = num_js_found + 1;
                Tij1j2(:,num_js_found) = Tijs(:,e);
            end
        end

        % finding real Tijtilde
        Tijtilde = zeros(nrs, low_deg);
        num_js_found = 0;
        for e = 1:size(edges)
            e_i = edges(e,1);
            if (e_i == ii)
                e_j = edges(e,2);
                num_js_found = num_js_found + 1;
                Tijtilde(:,num_js_found) = T(:, e_j) - T(:, e_i);
            end
        end

        [RitildeEst1,RitildeEst2,Qx_i]=recoverRitilde(R_i,Tijtilde);
        Qs_poc = [Qs_poc, Qx_i];
        disp('Again, one of the two residuals should be equal to zero')
        Ritilde = Qx_i' * [R_gt(:,:,ii); zeros(1,3)];
        disp(norm(RitildeEst1-Ritilde,'fro'))
        disp(norm(RitildeEst2-Ritilde,'fro'))
        disp('')

    end
end



