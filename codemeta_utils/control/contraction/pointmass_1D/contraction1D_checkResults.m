function [ ] = contraction1D_checkResults( kv, kd, beta, m )
%INPUT:
%kv - velocity gain > 0
%kd - position gain > 0
%beta - scalar > 0
%m - [3x ? ] matrix containing the m11, m12, m22 entries of the M matrix,
%   only the last row will be considered

%Check inputs
if (kv < 0 || kd < 0 || beta < 0)
    fprintf('ERROR: kv, kd, or beta is < 0\n');
    return;
end
if (size(m,1) ~= 3)
    fprintf('ERROR: m does not have 3 rows\n');
    return;
end

%Define standard system matrices
F=[0 1;-kd -kv];
m = m(:,end); %take the last entry, this should be the latest estimate
M = [m(1) m(2);m(2) m(3)];

%Check if the eigenvalues of LHS contraction is valid (all e-values < 0)
G=F'*M+M*F;
fprintf('eig(LHS) (should be < 0) = \n');
eig(G)
if (all(G <= -beta*M))
    fprintf('F''M+M*F <= -beta*M, TRUE\n');
else
    fprintf('F''M+M*F <= -beta*M, FALSE\n');
end

%Check if G is valid (all e-values < 0)
G=F'*M+M*F+beta*M;
fprintf('eig(G) (should be < 0) = \n');
eig(G)

end

