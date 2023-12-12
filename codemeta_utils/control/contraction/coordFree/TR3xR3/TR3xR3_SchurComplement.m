function [U,M_schur] = TR3xR3_SchurComplement(M)
% Compute the LDU decompsition from the Schur complment for the non-natural
% metric on TR^3xR^3.
% NOTE: On TR^3xR^3, the M matrix should be expanded to a [9x9] matrix using 
% the kronker product.
% INPUTS:
%   M := A symmetric [3x3] matrix where M(3,3) is the metric factor on R^3
%       and M(1:2,1:2) is the metric gain matrix on TR^3
% OUTPUTS:
%   U,M_schur := the resulting LDU decomposition where M = U'*M_schur*U

% Extract the block matrices of M
A = M(1:2,1:2);
B = M(1:2,3);
C = B';
D = M(3,3); % Since this is pos def., inv(D) exists

% Compute the L matrix
U = [eye(2) zeros(2,1);D\C eye(1)];
% Compute Diag matrix (V)
M_schur = [A-B*(D\C) zeros(2,1);zeros(1,2) D];
end

