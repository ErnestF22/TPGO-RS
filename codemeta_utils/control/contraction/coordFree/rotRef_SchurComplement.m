function [U,M_schur] = rotRef_SchurComplement(M)
% Compute the LDU decompsition from the Schur complment for the non-natural
% metric on TSO3xR3.
% NOTE: On TSO3xR3, the M matrix should be expanded to a [9x9] matrix using 
% the kronker product.
% INPUTS:
%   M := A symmetric [3x3] matrix where M(1,1) is the metric factor on SO(3)
%       and M(2:3,2:3) is the metric gain matrix on TSO(3)
% OUTPUTS:
%   U,M_schur := the resulting LDU decomposition where M = U'*M_schur*U

% Extract the block matrices of M
A = M(1,1);
B = M(1,2:3);
C = B';
D = M(2:3,2:3); % Since this is pos def., inv(D) exists

% Compute the L matrix
U = [1 zeros(1,2);D\C eye(2)];
% Compute Diag matrix (V)
M_schur = [A-B*(D\C) zeros(1,2);zeros(2,1) D];
end

