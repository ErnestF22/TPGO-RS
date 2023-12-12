function [ diff_Rk ] = rotBundle_RkDiff( Y, X, PU)
% Compute the differential of Rk (dK|Y) parameterize by (p,u) in TSO(3).
% The differential map is defined in Gudmudsson and Kappos, pg. 4 as:
% dRk|Y = d(exp_p o R_{-u} o tau)
% INPUTS:
%   Y := A point on the tangent bundle TSO(3) represented by a [6 x 3]
%       matrix. Y \in pi^(-1)(V). Y is assumed be in the form
%       [Rq; Rq*skew-symm];
%   X := A tangent vector in T_{p}_T_{u}SO(3) represented by a [6 x 3] 
%       matrix in the form [Rq*skew-symm; Rq*skew-symm]
%   PU := A point on the tangent bundle TSO(3) containing the point of
%       evaluation (p,u). NOTE: u is assumed to be in the form 
%       R*skew-symmetric matrix. NOTE: PU == Y(0) for connection map
% OUTPUTS:
%   diff_Rk := the differential of Rk rep. by a [3 x 3] matrix

% dRk|Y is a chain of the differential of each of its components


%% Differential of tau (sigma function in our case)
% Split Y into its components: Rq on SO(3) and Rq_vHat on T_q SO(3)
[Rq, Rq_vHat] = rotBundle_split(Y);
% Split Y_dot into its components: Rq_wHat on T_q SO(3) and Rq_zHat on
% T_(q,v)TSO(3). NOTE: Rq_wHat == Rq_vHat
[~, Rq_zHat] = rotBundle_split(X);
% Split PU into its componets: Rp on SO(3) and Rp_uHat on T_p SO(3)
[Rp, Rp_uHat] = rotBundle_split(PU);

dsigma = rot_parallelDiff(Rq, Rp, Rq_zHat);

%% Compute R_{-u} (same implmentation as rotBundle_Rk.m)
% the transported tangent vector: T_q SO(3) to T_p SO(3)
Rp_vHat = rot_parallel(Rq, Rp, Rq_vHat, 'toRotation');
% R_u is the translation of the vector Rp_vHat by vector Rp_uHat
R_u = Rp_vHat - Rp_uHat;

%% Differential of R_{-u}
dR_u = rot_rMinusUDiff(Rp, dsigma);

%% Differential of Exponential map and resulting diff_Rk
diff_Rk = rot3_expDiff(Rp, R_u, dR_u, 'tangent', 'method', 'closedForm');

end

