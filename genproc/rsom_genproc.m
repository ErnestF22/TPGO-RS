function [transf_out, rs_recovery_success] = rsom_genproc(T_gf, Tijs, edges, params, transf_initguess)
%RSOM_RS Rsom Manopt pipeline, with the addition of the Riemannian
%Staircase ("RS")

if ~exist('thresh','var')
  thr=1e-5;
end

rs_recovery_success = boolean(1);

nrs = size(T_gf, 1);
d = size(Tijs, 1);
N = size(T_gf, 2);
problem.sz = [nrs, d, N];

r0 = d+1;
problem_data.Tijs = Tijs;
problem_data.edges = edges;

tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
M = productmanifold(tuple);

% Setup the problem structure with manifold M and cost+grad functions.
problem.M = M;
problem_data.sz = [nrs, d, N];
problem.cost = @(x) cost_genproc(x, problem_data);
problem.grad = @(x) grad_genproc(x, problem_data);
problem.hess = @(x, u) hess_genproc(x, u, problem_data);

% checkgradient(problem);
% checkhessian(problem);

X = trustregions(problem);
T = X.T;
R = X.R;

cost_last = rsom_cost_base(X, problem_data);

for staircase_step_idx = r0:d*N+1
    problem_struct_next.sz = [staircase_step_idx, d, N];
    problem_struct_next.Tijs = Tijs;
    problem_struct_next.edges = edges;

    [Y_star, lambda, v] = rsom_pim_hessian_genproc( ...
        X, problem_struct_next, thr);
    disp("v") %just to remove unused variable warning
    disp(v)

    if lambda > 0
        disp("R, T eigenvals > 0: exiting staircase")
        break;
    end
    
    % next optimization iteration
    tuple_next.R = stiefelfactory(staircase_step_idx, d, N);
    tuple_next.T = euclideanfactory(staircase_step_idx, N);
    M_next = productmanifold(tuple_next);
    problem_next.M = M_next;    
    problem_next.cost = @(x) cost_genproc(x, problem_data); %!! problem_data is the same
    problem_next.grad = @(x) grad_genproc(x, problem_data);
    problem_next.hess = @(x, u) hess_genproc(x, u, problem_data);

    
    X = trustregions(problem_next, Y_star);
    T = X.T;
    R = X.R;

    disp("cost_last")
    disp(cost_last)
    cost_last = rsom_cost_base(X, problem_struct_next); 
    disp("cost_new")
    disp(cost_last)
    
    if rank(matStackH(Y_star.R))<staircase_step_idx
        break;
    end

end

if staircase_step_idx > d+1

    low_deg = 2; %TODO: not necessarily in more complex graph cases
    nodes_high_deg = params.node_degrees > low_deg;
    
    T_edges = make_T_edges(T, edges);
    
    RT_stacked_high_deg = [matStackH(R(:,:,nodes_high_deg)), T_edges];
    Qx = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);
    % RT_stacked_high_deg_poc = Qx * RT_stacked_high_deg;
    
    R_tilde = multiprod(repmat(Qx, 1, 1, sum(nodes_high_deg)), R(:,:,nodes_high_deg));
    
    R_out = zeros(d,d,N);
    R_out(:,:,nodes_high_deg) = R_tilde(1:d,:,:);
    
    
    % nodes_low_deg = ~nodes_high_deg;
    for node_id = 1:length(params.node_degrees)
        node_deg = params.node_degrees(node_id);
        if node_deg == low_deg
            fprintf("Running recoverRitilde() on node %g\n", node_id);
            R_i = R(:,:,node_id);
            
            check_Rb_params.node_degrees = params.node_degrees;
            check_Rb_params.nrs = staircase_step_idx - 1;
            check_Rb_params.d = d;

            [Tij1j2, Tij1j2_tilde] = make_Tij1j2s(node_id, R,T,Tijs,edges,check_Rb_params);

%             %computing Tij12
%             Tij1j2 = zeros(d,low_deg);
%             num_js_found = 0;
%             for e = 1:size(edges)
%                 e_i = edges(e,1);
%                 %         e_j = edges(e,2);
%                 if (e_i == node_deg)
%                     num_js_found = num_js_found + 1;
%                     Tij1j2(:,num_js_found) = Tijs(:,e);
%                 end
%             end
%     
%             % finding real Tijtilde
%             Tij1j2_tilde = zeros(staircase_step_idx-1, low_deg);
%             num_js_found = 0;
%             for e = 1:size(edges)
%                 e_i = edges(e,1);
%                 if (e_i == node_deg)
%                     e_j = edges(e,2);
%                     num_js_found = num_js_found + 1;
%                     Tij1j2_tilde(:,num_js_found) = T(:, e_j) - T(:, e_i);
%                 end
%             end
%             
%             check_Rb_params.nrs = staircase_step_idx - 1;
%             check_Rb_params.node_deg = node_deg;
    
            [RitildeEst1,RitildeEst2,Qx_i,Rb_i]=recoverRitilde(R_i,Tij1j2_tilde);
            check_Rb_ambiguity(node_id, Qx_i, Rb_i, R_i, Tij1j2, Tij1j2_tilde, check_Rb_params);
    %         Qs_poc = [Qs_poc, Qx_i];
    %         disp('Again, one of the two residuals should be equal to zero')
    %         Ritilde = Qx_i' * [R_gt(:,:,ii); zeros(1,3)];
    %         disp(norm(RitildeEst1-Ritilde,'fro'))
    %         disp(norm(RitildeEst2-Ritilde,'fro'))
            disp('')
            % TODO: how to decide between RitildeEst1,RitildeEst2??
            R_out(:,:,node_id) = RitildeEst1(1:d,:);
%             R_out = R_out;
            T_diffs = Qx * T_edges;
            T_out = edge_diffs_2_T(T_diffs(1:d,:), edges, N);
        end
    end
else
    R_out = R;
    T_out = T;
end

%checking that cost has not changed during "recovery"
X_out.T = T_out;
X_out.R = R_out;
cost_out = rsom_cost_base(X_out, problem_struct_next); 
disp("cost_out")
disp(cost_out)

transf_out = RT2G(R_out, T_out); %??

end %file function

%%
function Qx=align2d(v)
Q=fliplr(orthComplement(v));
Qx=flipud(orthCompleteBasis(Q)');
end

function RbEst=procrustesRb(c,q)
[U,~,V]=svd(c*q');
RbEst=U*diag([1 det(U*V')])*V';
end

function [RitildeEst1,RitildeEst2,Qx,RbEst]=recoverRitilde(Ritilde2,Tijtilde)
Qx=align2d(Tijtilde);
QxRitilde2Bot=Qx(3:4,:)*Ritilde2;
[U,~,~]=svd(QxRitilde2Bot,'econ');
c=U(:,2);

QLastRight=Qx(3:4,4)';

RbEst=procrustesRb(c,QLastRight');
RitildeEst1=Qx'*blkdiag(eye(2),-RbEst')*Qx*Ritilde2;
RitildeEst2=Qx'*blkdiag(eye(2),RbEst')*Qx*Ritilde2;
end

%% old methods
%TODO: Improve reprojection on SE(d)^N

% METHOD 1): Simply extract upper left submats
% transf_out = RT2G(R(1:d, 1:d, 1:N), T(1:d, 1:N));

% METHOD 2): SE-SYNC round solution (only for rotations!)
% R_out = matUnstackH( ...
%     round_solution_se_sync(matStackH(R), problem_struct_next));
% T_out = T(1:d, 1:N); %what to do??
% transf_out = RT2G(R_out, T_out);

% METHOD 3): CODEMETA lowRankLocalization_solution_extractProjection
% if (size(R, 1) > d)
%     [R_out, T_out] = lowRankLocalization_solution_extractProjection( ...
%     matStack(multitransp(R)) * T);
%     transf_out = RT2G(R_out, T_out);
% else
%     transf_out = RT2G(R, T);
% end

% METHOD 4 rots only): POC Rots only
% x_rs = matStackH(R);
% if staircase_step_idx > d+1
%     Q = POCRotateToMinimizeLastEntries(x_rs);
%     disp('x_rs=')
%     disp(x_rs)
%     disp('R=')
%     disp(R)
%     disp('x=Q*x_rs=')
%     disp(Q*x_rs)
%     x = Q*x_rs;
% else
%     x = x_rs;
% end
% transf_out = RT2G(matUnstackH(x(1:d, :)), T(1:d, 1:N));

% METHOD 4): POC with edges differences
% T_edges = make_T_edges(T, edges);
% x_rs = [matStackH(R), T_edges];
% if staircase_step_idx > d+1
%     Q = POCRotateToMinimizeLastEntries(x_rs);
%     disp('x_rs=')
%     disp(x_rs)
%     disp('R=')
%     disp(R)
%     disp('x=Q*x_rs=')
%     disp(Q*x_rs)
%     x = Q*x_rs;
%     if max(abs(x(d+1:end, :)), [], "all") > 1e-5
%         disp("max(abs(x(d+1:end, :)), [], ""all"") > 1e-5!!! " + ...
%             "-> x on SE(d)^N recovery failed")
%         rs_recovery_success = boolean(0);
%     end 
%     x_out = matUnstackH(x(1:d, 1:N*d));
%     T_diffs = x_rs(1:d, N*d+1:end);
%     [T_out, booleans_T] = edge_diffs_2_T(T_diffs, edges, N);
%     if min(booleans_T) < 1
%         disp("min(booleans_T) < 1!!! -> T recovery failed")
%         rs_recovery_success = boolean(0);
%     end 
% else
%     x_out = R;
%     T_out = T;
% end
% transf_out = RT2G(x_out, T_out); %??
