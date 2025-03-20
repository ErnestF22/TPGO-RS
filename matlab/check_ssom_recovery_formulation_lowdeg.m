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
Qy = POCRotateToMinimizeLastEntries(Y_stack);

R_tilde2_edges = multiprod(repmat(Qy, 1, 1, sum(nodes_high_deg)), R(:,:,nodes_high_deg));
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
            T_diffs_shifted = Qy * T_diffs; %this has last row to 0
            [a, b, Lambda_rec] = ssom_make_tij1j2s_edges(node_id, ...
                T_diffs_shifted, lambdas, tijs, edges, problem_data_local);

            %% (27), (28)
            disp("Qy * [R(:,:,node_id) * a *Lambda_rec, b]")
            disp([Qy * R(:,:,node_id) * a*Lambda_rec, b])
            disp("Checking (27), (28) -> diff between lhs and rhs should be 0")
            disp("max(abs(Qy * R(:,:,node_id) * a*Lambda_rec - b), [], ""all"")")
            disp(max(abs(Qy * Ri_tilde2 * a*Lambda_rec - b), [], "all"))

            %% (29)
            Qx = POCRotateToMinimizeLastEntries(b);
            cc = Qx * b;
            disp("cc")
            disp(cc)

            disp("Checking (29) -> last (p-d) rows should be 0")
            disp("max(abs(cc(d:end, :)), [], ""all"")")
            disp(max(abs(cc(d:end, :)), [], "all"))

            
            [RitildeEst1,RitildeEst2,Qx,Rb] = ...
                ssom_recoverRitilde(Qy * Ri_tilde2,b);

            %% (30)

            lhs_30 = Qx * b;

            disp("Checking (30) -> last (p-2) rows should be 0")
            disp("max(abs(lhs_30(d:end, :)), [], ""all"")")
            disp(max(abs(lhs_30(d:end, :)), [], "all"))

            %% (31a) -> RitildeEst1
            dd = Qx * RitildeEst1 * a * Lambda_rec;

            disp("dd")
            disp(dd)

            disp("Checking (31a) -> last (p-2) rows should be 0")

            disp("max(abs(dd(d:end, :)), [], ""all"")")
            disp(max(abs(dd(d:end, :)), [], "all"))

            %% (31b) -> RitildeEst2
            dd = Qx * RitildeEst2 * a * Lambda_rec;

            disp("dd")
            disp(dd)

            disp("Checking (31b) -> last (p-2) rows should be 0")

            disp("max(abs(dd(d:end, :)), [], ""all"")")
            disp(max(abs(dd(d:end, :)), [], "all"))

            %% (32)
            Qb = blkdiag(eye(low_deg), Rb);

            lhs32 = Qb * dd;

            disp("lhs32")
            disp(lhs32)

            disp("Checking (32) -> last (p-2) rows should be 0")
            
            disp("max(abs(lhs32(d:end, :)), [], ""all"")")
            disp(max(abs(lhs32(d:end, :)), [], "all"))

            %% (33a)
            lhs_33 = Qb * Qx * RitildeEst1 * a * Lambda_rec;
            rhs_33 = Qx * b;

            disp("Checking (33) -> lhs and rhs should be equal")

            disp("max(abs(lhs_33 - rhs_33)), [], ""all"")")
            disp(max(abs(lhs_33 - rhs_33), [], "all"))

            %% (33b)
            lhs_33 = Qb * Qx * RitildeEst2 * a * Lambda_rec;
            rhs_33 = Qx * b;

            disp("Checking (33) -> lhs and rhs should be equal")

            disp("max(abs(lhs_33 - rhs_33)), [], ""all"")")
            disp(max(abs(lhs_33 - rhs_33), [], "all"))
           

            %% (34a)
            lhs_34 = Qx' * Qb * Qx * RitildeEst1 * a * Lambda_rec;
            rhs_34 = b;

            disp("Checking (34a) -> lhs and rhs should be equal")

            disp("max(abs(lhs_34 - rhs_34)), [], ""all"")")
            disp(max(abs(lhs_34 - rhs_34), [], "all"))

            %% (34b)
            lhs_34 = Qx' * Qb * Qx * RitildeEst2 * a * Lambda_rec;
            rhs_34 = b;

            disp("Checking (34b) -> lhs and rhs should be equal")

            disp("max(abs(lhs_34 - rhs_34)), [], ""all"")")
            disp(max(abs(lhs_34 - rhs_34), [], "all"))

            %% (35a)
            lhs_35 = RitildeEst1 * a;
            
            disp("Checking (35a) -> last (p-3) rows should be 0")
            
            disp("max(abs(lhs_35(d+1:end, :)), [], ""all"")")
            disp(max(abs(lhs_35(d+1:end, :)), [], "all"))

            %% (36)

            %about the span of Q_x^{top}\transpose

            %% (37)
            lhs_37 = Qx' * Qb * Qx * RitildeEst1 * a * Lambda_rec;
            rhs_37 = b;

            disp("Checking (37) -> lhs and rhs should be equal")

            disp("max(abs(lhs_37 - rhs_37)), [], ""all"")")
            disp(max(abs(lhs_37 - rhs_37), [], "all"))
            

            %% (37) underbrace

            disp("Checking UNDERBRACE OF (37) -> lhs and rhs should be equal")

            disp("max(abs(Qx' * Qb * Qx * RitildeEst1 - Qy * Ri_tilde2)), [], ""all"")")
            disp(max(abs(Qx' * inv(blkdiag(eye(2),-Rb')) * Qx * RitildeEst1 - Qy * Ri_tilde2), [], "all"))
            disp(max(abs(Qx' * inv(blkdiag(eye(2),Rb')) * Qx * RitildeEst2 - Qy * Ri_tilde2), [], "all"))

            disp("")

            %% (38)

            % Rb = ssom_procrustesRb()

            % Qx' * blkdiag(eye(2), Rb') * Qx * Ri_tilde2

            %% (39)

            %direct consequence of (38)            


            %% (40)
            Qxtransp = Qx';
            Qtop = Qxtransp(1:d, :);
            Qbot = Qxtransp(d+1:end, :);
            Qbotleft = Qbot(:, 1:low_deg);
            Qbotright = Qbot(:, low_deg + 1:end);

            mathcalR = Qx * Qy * Ri_tilde2;
            Rtop = mathcalR(1:low_deg, :);
            Rbot = mathcalR(d:end, :);

            supposedly_zero_without_lemma = ...
                Qbotleft * Rtop + Qbotright * Rb' * Rbot; %supposedly null 

            disp("max(abs(supposedly_zero_without_lemma), [], ""all"")")
            disp(max(abs(supposedly_zero_without_lemma), [], "all"))
    
            %% (41)
            supposedly_zero_applying_lemma = Qbotright * Rb' * Rbot; %supposedly null

            disp("max(abs(supposedly_zero_applying_lemma), [], ""all"")")
            disp(max(abs(supposedly_zero_applying_lemma), [], "all"))

           
            

            


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