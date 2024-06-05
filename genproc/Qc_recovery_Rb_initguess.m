function [R_real_out, T_real_out, nodes_with_low_deg] = ...
    Qc_recovery_Rb_initguess(sz, edges, R, T, Tijs, node_degrees)
%QC_RECOVERY_RB_INITGUESS Summary of this function goes here
%   Detailed explanation goes here

nrs = sz(1);
d = sz(2);
N = sz(3);
p = nrs;

%% standard recovery method
T_edges = make_T_edges(T, edges);
x_rs = [matStackH(R), T_edges];
R_deg3 = matStackH(R(:,:,node_degrees == N-2));
x_R = matStackH(R_deg3);

Qa = POCRotateToMinimizeLastEntries(x_R);
disp('x_rs=')
disp(x_rs)
disp('R=')
disp(R)
disp('x=Q*x_R=')
disp(Qa*x_R)
x = Qa*x_R;
if max(abs(x(d+1:end, :)), [], "all") > 1e-5
    disp("max(abs(x(d+1:end, :)), [], ""all"") > 1e-5!!! " + ...
        "-> x on SE(d)^N recovery failed")
    rs_recovery_success = boolean(0);
end
% x_out = matUnstackH(x(1:d, 1:N*d));
% T_diffs = x_rs(1:d, N*d+1:end);
% [T_out, booleans_T] = edge_diffs_2_T(T_diffs, edges, N);
% if min(booleans_T) < 1
%     disp("min(booleans_T) < 1!!! -> T recovery failed")
%     rs_recovery_success = boolean(0);
% end

% transf_out = RT2G(x_out, T_out); %??


%% find Qb

%set up data

%trying with i=1
%its edges are j1 = 2, j2 = 3


low_deg = 2;
nodes_with_low_deg = node_degrees==2;

Qas = zeros(p,p,N);

for ii = 1:N
    if ~nodes_with_low_deg(ii)
        continue;
    end
    j1 = 2;
    j2 = 3;
    R_i = R(:,:,1);
    idx_i_j1 = find(ismember(edges, [ii, j1], "rows"));
    idx_i_j2 = find(ismember(edges, [ii, j2], "rows"));
    T_i_j1_j2 = [Tijs(:,idx_i_j1), Tijs(:,idx_i_j2)];
    T_i_j1_j2_tilde = - [T(:,ii) - T(:,j1), T(:,ii) - T(:,j2)]; % !! -

    disp("max(abs(R_i * T_i_j1_j2 - T_i_j1_j2_tilde),[], ""all"")")
    disp(max(abs(R_i * T_i_j1_j2 - T_i_j1_j2_tilde),[], "all"))
    disp("max(abs(Qa * R_i * T_i_j1_j2 - Qa* T_i_j1_j2_tilde),[], ""all"")")
    disp(max(abs(Qa * R_i * T_i_j1_j2 - Qa* T_i_j1_j2_tilde),[], "all"))

    % Create the problem structure.
    manifold = stiefelfactory(p-2,p-2);
    problem_Qbnn.M = manifold;

    % Define the problem cost function and its Euclidean gradient.
    problem_Qbnn.cost  = @(x) mycost_Qbnn(Qa, x, R_i, T_i_j1_j2, T_i_j1_j2_tilde);
    problem_Qbnn = manoptAD(problem_Qbnn);

    % Numerically check gradient consistency (optional).
    % checkgradient(problem_Qbnn);

    % Solve.
    options.maxiter = 100;
    [Rb, xcost, info, options] = trustregions(problem_Qbnn,[],options);

    Qb = blkdiag(eye(low_deg), Rb);
    disp("Qb")
    disp(Qb)

    disp("max(abs(Qb * Qa * R_i * T_i_j1_j2 - Qa* T_i_j1_j2_tilde),[], ""all"")")
    disp(max(abs(Qb * Qa * R_i * T_i_j1_j2 - Qa* T_i_j1_j2_tilde),[], "all"))

    Qa_i = POCRotateToMinimizeLastEntries(x_R);
    Qas(:,:,ii) = Qa_i;
    disp("max(abs(Qa' * Qb * Qa * R_i * T_i_j1_j2 - T_i_j1_j2_tilde), [], ""all"")")
    disp(max(abs(Qa' * Qb * Qa * R_i * T_i_j1_j2 - T_i_j1_j2_tilde), [], "all"))
end
%% check on all Rs

% step1 = x; % x = Qa*x_R
% disp("max(abs(step1(end, :)), [], ""all"")")
% disp(max(abs(step1(end, :)), [], "all"))
%
% step2 = Qb * step1;
% disp("max(abs(step2(end, :)), [], ""all"")")
% disp(max(abs(step2(end, :)), [], "all"))
%
% disp("max(abs(step1-step2), [], ""all"")")
% disp(max(abs(step1-step2), [], "all"))

% step3 = Qa' * Qb * step2;
% disp("max(abs(step3(end, :)), [], ""all"")")
% disp(max(abs(step3(end, :)), [], "all"))

% x_R = matStackH(R);
% disp("Qa * x_R")
% disp(Qa * x_R)
%
% disp("Qb * Qa * x_R")
% disp(Qb * Qa * x_R)
%
% disp("Qa' * Qb * Qa * x_R")
% disp(Qa' * Qb * Qa * x_R)

%% actual Qas
disp("actual Qas")
for ii = 1:N
    if ~nodes_with_low_deg(ii)
        continue;
    end

    Qa_i = POCRotateToMinimizeLastEntries(R(:,:,ii));
    fprintf("actual Qa_i for i = %g -> [Qa_i, Qa_i*R(:,:,ii)]\n", ii)
    disp([Qa_i, Qa_i*R(:,:,ii)])
end

%% more complex recovery method (Qa, Qb, Qcd, Qcdd, linsyst+procr initguess)

rots_recovered = zeros(p,d,N);
for ii = 1:N
    if ~nodes_with_low_deg(ii)
        continue;
    end
    %     tuple.qc = stiefelfactory(p,p);
    %     tuple.rb = stiefelfactory(p-low_deg,p-low_deg);
    tuple.qc = rotationsfactory(p);
    tuple.rb = rotationsfactory(p-low_deg);
    problem_qcrb.M = productmanifold(tuple);

    fprintf("Now running recovery procedure on node %g\n", ii);
    Ri = R(:,:,ii);
    Qa_i = POCRotateToMinimizeLastEntries(Ri);

    [Rb_initguess] = find_Rb_initguess(low_deg,p,d,Qa_i,R_i);
    disp('Rb_initguess')
    disp(Rb_initguess)

    Qcd_i = eye(p); % Hp) to be checked!

    % Define the problem cost function and its Euclidean gradient.
    problem_qcrb.cost = @(x) mycost_qcrb( ...
        x, low_deg, Qa_i, Qcd_i, Ri, T_i_j1_j2, T_i_j1_j2_tilde);
    problem_qcrb.egrad = @(x) myegrad_qcrb( ...
        x, low_deg, Qa_i, Qcd_i, Ri, T_i_j1_j2, T_i_j1_j2_tilde);

    % Numerically check gradient consistency (optional)
    %     x_chkgrad.qc = eye(4);
    %     x_chkgrad.rb = make_rand_stiefel_3d_array(2,2,1);
    %     checkgradient(problem_qcrb, x_chkgrad);
    checkgradient(problem_qcrb);

    % Solve providing initguess.
    initguess_j.qc = make_rand_stiefel_3d_array(p,p,1);
    if det(initguess_j.qc) < 0
        initguess_j.qc(:,1) = - initguess_j.qc(:,1);
    end
    initguess_j.rb = Rb_initguess;
    if det(initguess_j.rb) < 0
        initguess_j.rb(:,1) = - initguess_j.rb(:,1);
    end

    disp('initguess_qc -> check_is_rotation()')
    disp(check_is_rotation(initguess_j.qc))
    disp('initguess_rb -> check_is_rotation()')
    disp(check_is_rotation(initguess_j.rb))

    %     problem_qcrb = manoptAD(problem_qcrb);

    [qcrb_out_i, xcost, info, options] = ...
        trustregions(problem_qcrb,initguess_j,options);

    % Solve WITHOUT initguess.
    %     [qcrb_out_i, xcost, info, options] = ...
    %         trustregions(problem_qcrb,[],options);

    rot_i_recov = qcrb_out_i.qc' * Qcd_i * Ri;
    disp("rot_i_recov")
    disp(rot_i_recov)
    rots_recovered(:,:,ii) = rot_i_recov;
end

%check if recovery was actually a success
for ii = 1:N
    if ~nodes_with_low_deg(ii)
        continue;
    end
    fprintf('check_is_rotation(R_final(:,:,%g))\n',ii)
    disp(check_is_rotation(rots_recovered(1:d,:,ii)))
end


R_real_out = rots_recovered;
T_real_out = zeros(d,N); %TODO: Fix this!!

end %file function


%% cost functions

% function c_out = mycost_Qbnn(Qa, Rb, Ri, Ti, Ti_tilde)
%     p = size(Qa, 1);
%     Qb = eye(p);
%     Qb(p-1:end, p-1:end) = Rb;
%     c = Qa' * Qb * Qa * Ri;
%     c_out = norm(c(end,:));
% end

function c_out = mycost_Qbnn(Qa, Rb, Ri, Ti, Ti_tilde)
p = size(Qa, 1);
Qb = eye(p);
Qb(p-1:end, p-1:end) = Rb;
c = Qa' * Qb * Qa * Ri * Ti - Ti_tilde;
c_out = norm(c(:));
end

