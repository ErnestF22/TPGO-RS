function [U,M_schur] = rotBundle_SchurComplement(M)
% Compute the LDU decompsition from the Schur complment for the non-natural
% metric on TSO3.
% NOTE: On TSO3, the M matrix should be expanded to a [6x6] matrix using 
% the kronker product.
% INPUTS:
%   M := A symmetric [2x2] matrix representing the nonnatural metric tensor
% OUTPUTS:
%   U,M_schur := the resulting LDU decomposition where M = U'*M_schur*U

A=M(1,1); B = M(1,2);
C = M(2,1); D = M(2,2);
M_schur = [A-B*inv(D)*C 0;0 D];
U = [1 0; inv(D)*C 1];
end

