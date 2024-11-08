% Supplemetary material for the SIMAX manuscript
%
% "A matrix-algebraic algorithm for the Riemannian logarithm on the 
%    Stiefel manifold under the canonical metric"

function [U1] = Stiefel_Exp_supp(U0, Delta)
%-------------------------------------------------------------
%@author: Ralf Zimmermann, IMADA, SDU Odense
% Input arguments      
%   U0    : base point on St(n,p)
%   Delta : tangent vector in T_U0 St(n,p)
% Output arguments
%   U1    : Exp^{St}_U0(Delta), 
%-------------------------------------------------------------
% get dimensions
[n,p] = size(U0);
A = U0'*Delta;                          % horizontal component
K = Delta-U0*A;                             % normal component
[Qe,Re] = qr(K, 0);                   % qr of normal component
%OBS. qr(A,0) is equivalent to setting both "econ" and "vector" options.
% matrix exponential
MNe = expm([[A, -Re'];[Re, zeros(p)]]);
U1 = [U0, Qe]*MNe(:,1:p);
return;
end