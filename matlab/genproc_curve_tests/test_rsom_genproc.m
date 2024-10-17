function problem=test_rsom_genproc()

testdata = testNetwork_som(3); %4 is the default option

nrs = 3;
d = 3;
N = testdata.NNodes;
edges = (testdata.E);
% num_edges = testdata.NEdges;

sz=[nrs,d,N];

Tijs = G2T(testdata.gijtruth);
R_gf = G2R(testdata.gitruth);
T_gf = G2T(testdata.gitruth);

T_gf_stief = cat_zero_row(T_gf);

[P, frct] = make_step1_p_fct(T_gf_stief, Tijs, edges);
[LR, PR, BR] = make_LR_PR_BR_noloops(R_gf, Tijs, edges);


problem=struct("R", R_gf, "T", T_gf,  ...
    "Tijs", Tijs, "edges", testdata.E,...
    "sz", sz, ...
    'P', P, 'frct', frct, ...
    'PR', PR, 'LR', LR, 'BR', BR);



problem.cost=@(x) rsom_cost_rot_stiefel(x,problem);
problem.egrad=@(x) rsom_egrad_rot_stiefel(x,problem);
problem.grad=@(x) rsom_rgrad_rot_stiefel(x,problem);
problem.ehess=@(x,u) rsom_ehess_rot_stiefel(x,u,problem);
problem.hess=@(x,u) rsom_rhess_rot_stiefel(x,u,problem);


% problem.gs1 = @(x) grad_step1(x, problem);
% problem.hrt = @(x) hrt(x, problem);

problem.gs2 = @(x) grad_step2(x, problem);
problem.htr = @(x) htr(x, problem);


end %file function


function gs1 = grad_step1(x, problem)
    d = problem.sz(2);
    Tijs = problem.Tijs;
    R = problem.R;
    % P
    P = zeros(nrs, d*N);
    idx_col_p = reshape(1:d*N, [], N)';
    num_edges = size(edges,1);
    for e = 1:num_edges
        ii = edges(e,1);
        jj = edges(e,2); 
    %     R_i = R_gf(:,:,ii);
        T_j = x(:, jj);
        T_i = x(:, ii);
        T_ij = Tijs(:,e);
        P_e = 2 * (T_i * T_ij' - T_j * T_ij');
        P(:, idx_col_p(ii, :)) = ...
            P(:, idx_col_p(ii, :)) + P_e;
    end
    eg=matUnstackH(P,d);
    stiefel_tangentProj(R,eg);

    


end

function hrt_out = hrt(x, problem)
    R = problem.R;
    Tijs = problem.Tijs;

    W = zeros(size(R));
    num_edges = size(edges, 1);
    for e = 1:num_edges
        ii = edges(e,1);
        jj = edges(e,2); 
        v_ti = x(:,ii);
        v_tj = x(:,jj);
        T_ij = Tijs(:, e);
        w_ij = 2 * (v_ti - v_tj)*T_ij';
        W(:,:,ii) = W(:,:,ii) + w_ij;
    end
    hrt_out = stiefel_tangentProj(R, W);

end

function gs2 = grad_step2(x, problem)
%    gs2 = rsom_rgrad_transl_stiefel(x, problem);
    nrs = size(problem.R,1);
    N = size(problem.R, 3);
    num_edges = size(problem.edges, 1);
    %
    LR = zeros(N,N);
    PR = zeros(N,nrs);
%     BR_const = zeros(d,d);
    %     
    for ee = 1:num_edges
        BIJ = zeros(N,1);
        ii = problem.edges(ee, 1);
        jj = problem.edges(ee, 2);
        BIJ(ii) = 1;
        BIJ(jj) = -1;
        Tij = problem.Tijs(:, ee);
        Ri = problem.R(:,:,ii);
        %LR
        LR = LR + BIJ * BIJ';
        %PR
        PR = PR + 2 * BIJ * Tij' * Ri';
        %BR
%         BR_const = BR_const + Tij * Tij';
    end

    T = problem.T;
    gs2 = T*(problem.LR+problem.LR')+(problem.PR)';
end

function htr_out = htr(x, problem)
    %code from compute_htr.m
    %
    htr_out = zeros(problem.sz(1), problem.sz(3));
    N = problem.sz(3);
    num_edges = size(problem.edges, 1);
    %
    for e = 1:num_edges
        ii = problem.edges(e,1);
        jj = problem.edges(e,2); 
        v_ri = problem.R(:,:,ii);
%         v_rj = u.R(:,:,jj);
        %
        BIJ = zeros(N,1);
        BIJ(ii) = 1;
        BIJ(jj) = -1;
        %
        T_ij = problem.Tijs(:, e);
        w_ij = 2 * BIJ * T_ij' * v_ri';
        htr_out = htr_out + w_ij';
    end
end

