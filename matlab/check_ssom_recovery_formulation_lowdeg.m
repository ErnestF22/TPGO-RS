function check_ssom_recovery_formulation_lowdeg

load('data/ssom_recovery/ws2.mat', 'X_manopt_out')
load('data/ssom_recovery/ws2.mat', 'problem_data')
load('data/ssom_recovery/ws2.mat', 'problem_data_next')
load('data/ssom_recovery/ws2.mat', 'N')

R = X_manopt_out.R;
T = X_manopt_out.T;
lambdas = X_manopt_out.lambda;
edges = problem_data.E;
tijs = problem_data.tijs;

T_diffs = make_T_edges(T,edges);

node_degrees = sum(problem_data.A, 2);
low_deg = 2; %This is still an assumption for SSOM
nodes_high_deg = node_degrees > low_deg; % all nodes are high deg in this test
nodes_low_deg = ~nodes_high_deg;

%% from deg>=3 recovery procedure

% Y_stack = [matStackH(R), T_diffs];

% Qy = svd(Y_stack * Y_stack'); %last p-d rows are indeed 0
% Qx_edges = POCRotateToMinimizeLastEntries(Y_stack);
% R_tilde2_edges_high_deg = multiprod(repmat(Qx_edges, 1, 1, sum(nodes_high_deg)), R(:,:,nodes_high_deg));

% R(:,:,nodes_low_deg) = multiprod (repmat(Qx_edges, 1, 1, sum(nodes_low_deg)), R(:,:,nodes_low_deg));

%%

p = size(R,1);
d = size(R,2);

node_degrees = sum(problem_data.A, 2);
low_deg = 2; %This is still an assumption for SSOM
nodes_high_deg = node_degrees > low_deg; % all nodes are high deg in this test
nodes_low_deg = ~nodes_high_deg;




%% (24)
Y_stack = [matStackH(R(:,:,nodes_high_deg)), T_diffs];

%% (25)
Qy_25 = svd(Y_stack * Y_stack'); %last p-d rows are indeed 0
disp("Qy_25")
disp(Qy_25)

disp("max(abs(Qy(d+1:end, :)), [], ""all"")")
disp(max(abs(Qy_25(d+1:end, :)), [], "all"))

%% (26)
% T_diffs_shifted = Qy' * T_diffs


%%
Qx_edges = POCRotateToMinimizeLastEntries(Y_stack);

R_tilde2_edges = multiprod(repmat(Qx_edges, 1, 1, sum(nodes_high_deg)), R(:,:,nodes_high_deg));
disp("R_tilde2_edges")
disp(R_tilde2_edges)

disp("max(abs(R_tilde2_edges(d+1:end, :)), [], ""all"")")
disp(max(abs(R_tilde2_edges(d+1:end, :)), [], "all"))

R_recovered = zeros(d,d,N);
R_recovered(:,:,nodes_high_deg) = R_tilde2_edges(1:d,:,:);

T_recovered = zeros(d,d,N);

if ~any(nodes_low_deg)
    disp("See other test case for that!")
else
    for node_id = 1:length(node_degrees)
        node_deg = node_degrees(node_id);        

        if (node_id > 1)
            break; %TODO: remove this early break later
        end

        if node_deg == low_deg           

            fprintf("Running recoverRitilde() on node %g\n", node_id);
            Ri_tilde2 = R(:,:,node_id);

            cost_manopt_output = ssom_cost(X_manopt_out, problem_data_next); 
            disp("cost_manopt_output")
            disp(cost_manopt_output)
    

            problem_data_local = problem_data;
            problem_data_local.d = problem_data.sz(2);
            problem_data_local.node_degrees = node_degrees;
            T_diffs_shifted = Qx_edges * T_diffs; %this has last row to 0
            [a, b, Lambda_rec] = ssom_make_tij1j2s_edges(node_id, ...
                T_diffs_shifted, lambdas, tijs, edges, problem_data_local);

            %% (27), (28)
            disp("[R(:,:,node_id) * a *Lambda_rec, b]")
            disp([R(:,:,node_id) * a*Lambda_rec, b])
            disp("Checking (27), (28) -> residual should be 0")
            disp("max(abs(R(:,:,node_id) * a*Lambda_rec - b), [], ""all"")")
            disp(max(abs(Ri_tilde2 * a*Lambda_rec - b), [], "all"))

            %% (29)
            Qy = POCRotateToMinimizeLastEntries(b);
            cc = Qy * b;
            disp("cc")
            disp(cc)

            disp("Checking (29) -> residual should be 0")
            disp("max(abs(cc(d:end, :)), [], ""all"")")
            disp(max(abs(cc(d:end, :)), [], "all"))

            %% (30)
            dd = Qy * Ri_tilde2 * a * Lambda_rec;

            disp("dd")
            disp(dd)

            disp("Checking (30) -> residual should be 0")

            disp("max(abs(dd(d:end, :)), [], ""all"")")
            disp(max(abs(dd(d:end, :)), [], "all"))

            %% (31)
            Rb = make_rand_stiefel_3d_array(p-low_deg, p-low_deg, 1);
            Qb = blkdiag(eye(low_deg), Rb);

            lhs31 = Qb * dd;

            disp("lhs31")
            disp(lhs31)

            disp("Checking (31) -> residual should be 0")
            
            disp("max(abs(lhs31(d:end, :)), [], ""all"")")
            disp(max(abs(lhs31(d:end, :)), [], "all"))

            % %% (32)
            [RitildeEst1,RitildeEst2,Qx,Rb] = ...
                ssom_recoverRitilde(Qy * Ri_tilde2,b);
            lhs_32 = Qb * dd;
            rhs_32 = Qx * b;

            disp("Checking (32) -> residual should be 0")

            disp("max(abs(lhs32(d:end, :)), [], ""all"")")
            disp(max(abs(lhs_32 - rhs_32), [], "all"))
            % 
            % %% (33)
            % Qx = make_rand_stiefel_3d_array(p, p, 1);
            % lhs_32 = Qb * Qx * Ri_tilde2 * a;
            % rhs_32 = Qx * b;
            % 
            % disp("Checking (32) -> residual should be 0")
            % 
            % disp("max(abs(lhs31(d:end, :)), [], ""all"")")
            % disp(max(abs(lhs31(d:end, :)), [], "all"))
            
            
            Qxtransp = Qx';
            Qbot = Qxtransp(d+1:end, :);
            % mathcalR = Qy * Ri_tilde2;
            mathcalR = Qx * Ri_tilde2;
            Rbot = mathcalR(d:end, :);
            Qbotright = Qbot(:, d:end);

            Qbotleft = Qbot(:, 1:2); %supposedly leads to null term
            Rtop = mathcalR(1:2, :); %supposedly leads to null term

            %% (40)

            supposedly_zero_without_lemma = ...
                Qbotleft * Rtop + Qbotright * Rb' * Rbot; %supposedly null 

            disp("max(abs(supposedly_zero_without_lemma), [], ""all"")")
            disp(max(abs(supposedly_zero_without_lemma), [], "all"))
    
            %% (41)
            supposedly_zero_applying_lemma = Qbotright * Rb' * Rbot; %supposedly null

            disp("max(abs(supposedly_zero_applying_lemma), [], ""all"")")
            disp(max(abs(supposedly_zero_applying_lemma), [], "all"))

           
            %% (37)
            tmp = blkdiag(eye(2,2), Rb);
            lhs_37 = Qx' * tmp * Qx * RitildeEst1 * a * Lambda_rec;
            rhs_37 = b;

            disp("[lhs_37, rhs_37]")
            disp([lhs_37, rhs_37])

            %% (38)

            % Rb = ssom_procrustesRb()

            % Qx' * blkdiag(eye(2), Rb') * Qx * Ri_tilde2

            
            %% ssom_recoverRitilde()

            % [RitildeEst1,RitildeEst2,~,~] = ...
            %         ssom_recoverRitilde(Qx_edges * Ri_tilde2 , Qx * b);
            % 
            % 
            % disp("RitildeEst1")
            % disp(RitildeEst1)
            % disp("RitildeEst2")
            % disp(RitildeEst2)


        end
    end
end

% disp("R_recovered")
% disp(R_recovered)
% disp("T_recovered")
% disp(T_recovered)



end %file function

function [RitildeEst1, RitildeEst2, Qx, RbEst] = ssom_recoverRitilde(Ritilde2,Tijtilde)
Qx=align2d(Tijtilde);
QxRitilde2Bot=Qx(3:end,:)*Ritilde2;
[U,~,~]=svd(QxRitilde2Bot,'econ');
c=U(:,2);

QLastRight=Qx(3:end,4)';

RbEst=ssom_procrustesRb(c,QLastRight');
RitildeEst1=Qx'*blkdiag(eye(2),-RbEst')*Qx*Ritilde2;
RitildeEst2=Qx'*blkdiag(eye(2),RbEst')*Qx*Ritilde2;
% RitildeEst1 = rand(size(Ritilde2));
% RitildeEst2 = rand(size(Ritilde2));
end
    
function RbEst=ssom_procrustesRb(c,q)
[U,~,V]=svd(c*q');
diagmat = eye(size(U,1));
diagmat(end, end) = det(U*V');
RbEst=U*diagmat*V';
end