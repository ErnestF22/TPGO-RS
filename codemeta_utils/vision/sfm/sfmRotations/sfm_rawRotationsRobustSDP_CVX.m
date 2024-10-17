function X = sfm_rawRotationsRobustSDP_CVX(M,varargin)

% This is the implemention of the following papers
% L. Wang, A. Singer, ``Exact and Stable Recovery of Rotations for Robust Synchronizationâ€?
% J. Saunderson, et al, Semidefinite relaxations for optimization problems over rotation matrices
% which solves
% min \sum ||Mij-Xij||_2
% st: X is PSD, X_ii = I
% if proj2so3conv == true
% X_ij \in conv(SO(3))

proj2so3conv = false;
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'proj2so3conv'
            proj2so3conv = true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

A{1} = diag([1 1 -1 -1]);
A{2} = zeros(4,4); A{2}(1,4) = -1; A{2}(2,3) = 1; A{2} = A{2} + A{2}';
A{3} = zeros(4,4); A{3}(1,3) = 1; A{3}(2,4) = 1; A{3} = A{3} + A{3}';
A{4} = zeros(4,4); A{4}(1,4) = 1; A{4}(2,3) = 1; A{4} = A{4} + A{4}';
A{5} = diag([1 -1 1 -1]);
A{6} = zeros(4,4); A{6}(1,2) = -1; A{6}(3,4) = 1; A{6} = A{6} + A{6}';
A{7} = zeros(4,4); A{7}(1,3) = -1; A{7}(2,4) = 1; A{7} = A{7} + A{7}';
A{8} = zeros(4,4); A{8}(1,2) = 1; A{8}(3,4) = 1; A{8} = A{8} + A{8}';
A{9} = diag([1 -1 -1 1]);
B = zeros(16,9);
for i = 1:9
    B(:,i) = A{i}(:);
end

n = size(M,1)/3;

%warning off
cvx_begin sdp quiet
    variable X(3*n,3*n) symmetric
    maximize trace(M'*X)
    subject to
    X >= 0;
    for i = 1:n
        X(3*i-2:3*i,3*i-2:3*i) == eye(3)
        if proj2so3conv
            for j=i+1:n
                reshape(B*reshape(X(3*i-2:3*i,3*j-2:3*j),9,1),4,4) + eye(4) >= 0;
            end
        end
    end 
cvx_end

% for i = 1:n-1
%     for j = i+1:n
%         X(3*i-2:3*i,3*j-2:3*j) = projectR(X(3*i-2:3*i,3*j-2:3*j));
%         X(3*j-2:3*j,3*i-2:3*i) = X(3*i-2:3*i,3*j-2:3*j)';
%     end
% end

end
