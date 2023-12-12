function [ Rk ] = rotBundle_Rk( Y, PU )
% Given a Y in TSO(3), return the R_k mapping (the 'integral' of the
% connection mapping K(Z) used to define the horizontal subspace.)
% NOTE: the algorithm used is described in Gudmundsson and Kappos pg. 4
% INPUTS:
%   Y := A point on the tangent bundle TSO(3) represented by a [6 x 3]
%       matrix. Y \in pi^(-1)(V)
%   PU := A point on the tangent bundle TSO(3) containing the point of
%       evaluation (p,u). Note: u is assumed to be in the form 
%       R*skew-symmetric matrix. NOTE: PU == Y(0) for connection map
% OUTPUTS:
%   R_k := The mapping of Y from TSO(3) to SO(3)

%% Parallel Transport: sigma function (tau in paper)
% Parallel transport Y to PU
% Split Y into its components: Rq on SO(3) and Rq_vHat on T_q SO(3)
[Rq, Rq_vHat] = rotBundle_split(Y);
% Split PU into its componets: Rp on SO(3) and Rp_uHat on T_p SO(3)
[Rp, Rp_uHat] = rotBundle_split(PU);

% the transported tangent vector: T_q SO(3) to T_p SO(3)
Rp_vHat = rot_parallel(Rq, Rp, Rq_vHat, 'toRotation');

%% Translation: R_-u on T_p SO(3)
% R_u is the translation of the vector Rp_vHat by vector Rp_uHat
R_u = Rp_vHat - Rp_uHat;

%% Exponential map:
% maps R_u to SO(3) using the exponential map
Rk = rot_exp(Rp, R_u);

end

