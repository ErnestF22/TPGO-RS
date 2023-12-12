function [ diff_RminusU ] = rot_rMinusUDiff( R, vp_dot )
% Return the differential of the R mius mapping (dR_{-u}|Rp)
% INPUTS:
%   R := Point on SO(3) where the tangent vector u is defined rep. by 
%       [3 x 3] matrix
%   vp_dot := time derivative of the tangent vector that was parallel
%       transported from Rq to Rp, rep. by [3 x 3] matrix, assumed to in
%       the form of R*skew-symm. matrix
% OUTPUTS:
%   diff_RminusU := differential of the R_{-u} mapping rep. by [3 x 3]
%       matrix

% The mapping R_{-u} := vp_dot - Rp*uHat in T_p SO(3), where uHat is a
% skew-symm matrix, thus the differential (dR_{-u} = vp_dot.
% REMINDER: vp_dot = Rp*skew-symm matrix.

diff_RminusU = vp_dot;

end

