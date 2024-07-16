function test_POCdegree2Align_gt
% From work with Ernesto
% See equations (41)-(45) in our working document
flagRandomizeQx=true;
flagPerturbMeasurements=false;

resetRands()

% Generate relative translations (with the last row already set to zero)
N = 6;
load('Qbnn_data/testdata.mat','edges')
adj_mat = make_adj_mat_from_edges(edges, N);
testdata = testNetwork_adj(3, adj_mat, 'banded', 3);
G = graph(make_adj_mat_from_edges(testdata.E,N));
figure(1000)
plot(G)
title('graph')


load('Qbnn_data/R_gt.mat','R_gt')
% load('Qbnn_data/testdata.mat','T')
load('Qbnn_data/testdata.mat','Tijs')

% nrs = size(R, 1);
% p = nrs;

R = R_gt;

d = size(R, 2);

% setup node degrees
low_deg = 2;
% node_degrees = [2,2,4,3,3,4];
% nodes_high_deg = node_degrees == 3 | node_degrees == 4;
% nodes_low_deg = ~nodes_high_deg;

%% Running procedure on node ii = 1

ii = 1;


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

low_deg = 2;

for deg = node_degrees
    ii = deg;
    if ii == low_deg
        fprintf("Running recoverRitilde() on node %g\n", ii);
        R_i = R(:,:,ii);
        
        
        
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
        
        Ritilde = [R_gt(:,:,ii);zeros(1,3)];
        % Tijtilde=Ritilde*Tij1j2;
        
        % finding real Tijtilde
        Tijtilde = zeros(nrs, low_deg);
        num_js_found = 0;
        for e = 1:size(edges)
            e_i = edges(e,1);
            if (e_i == ii)
                e_j = edges(e,2);
                num_js_found = num_js_found + 1;
                Tijtilde(:,num_js_found) = T(:, e_j) - T(:, ii);
            end
        end
        
        
        Qx=align2d(Tijtilde);
        
        
        Rb=randrot_som(2);
        Qb=blkdiag(eye(2),Rb);        
        % Check that Qx'*Qb*Qx leaves Tijtilde invariant
        % And that the new Ritilde generated from the ambiguity
        % still gives the same measurements.
        Ritilde2=Qx'*Qb*Qx*Ritilde;
        
        [RitildeEst1,RitildeEst2]=recoverRitilde(Ritilde2,Tijtilde);
        disp('Again, one of the two residuals should be equal to zero')
        Ritilde = [R_gt(:,:,ii); zeros(1,3)];
        disp(norm(RitildeEst1-Ritilde,'fro'))
        disp(norm(RitildeEst2-Ritilde,'fro'))

    end
end



%%
function Qx=align2d(v)
Q=fliplr(orthComplement(v));
Qx=flipud(orthCompleteBasis(Q)');

function RbEst=procrustesRb(c,q)
[U,~,V]=svd(c*q');
RbEst=U*diag([1 det(U*V')])*V';

function [RitildeEst1,RitildeEst2]=recoverRitilde(Ritilde2,Tijtilde)
Qx=align2d(Tijtilde);
QxRitilde2Bot=Qx(3:4,:)*Ritilde2;
[U,~,~]=svd(QxRitilde2Bot,'econ');
c=U(:,2);

QLastRigh=Qx(3:4,4)';

RbEst=procrustesRb(c,QLastRigh');
RitildeEst1=Qx'*blkdiag(eye(2),-RbEst')*Qx*Ritilde2;
RitildeEst2=Qx'*blkdiag(eye(2),RbEst')*Qx*Ritilde2;
