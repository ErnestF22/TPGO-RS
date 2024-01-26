function [LR, PR, BR_const] = make_LR_PR_BR(R_gf, Tijs, edges)
%MAKE_LR_PR_BR Summary of this function goes here
%   Detailed explanation goes here

nrs = size(R_gf, 1);
d = size(R_gf, 2);
N = size(R_gf, 3);

num_edges = size(edges, 1);

LR = zeros(N,N);
PR = zeros(N,nrs);
BR_const = zeros(d,d);


for ee = 1:num_edges
    BIJ = zeros(N,1);
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    BIJ(ii) = 1;
    BIJ(jj) = -1;
    Tij = Tijs(:, ee);
    Ri = R_gf(:,:,ii);
    %LR
    LR = LR + BIJ * BIJ';
    %PR
    PR = PR + 2 * BIJ * Tij' * Ri';
    %BR
    BR_const = BR_const + Tij * Tij';
end


end %file function

