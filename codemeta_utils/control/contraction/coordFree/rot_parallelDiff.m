function [ diff_d ] = rot_parallelDiff( Rq, Rp, d_dot )
% Return the differential of the parallel transport (dT_Rp^Rq|Rq) from Rq
% to Rp.
% INPUTS:
%   Rq := Starting point on SO(3) rep. by [3 x 3] matrix
%   Rp := Ending point on SO(3) rep. by [3 x 3] matrix
%   d_dot := time derivative of the tangent vector in T_q SO(3), assumed to
%       be in the form or Rq*Skew-symmetric matrix, rep. by [3 x 3] matrix
% OUTPUTS:
%   diff_d := differential of the parallel transport rep. by [3 x 3] matrix

% The tangent vector is the only time-dependent variable, thus the
% differential is exactly the parallel transport mapping of d_dot

diff_d = rot_parallel(Rq, Rp, d_dot, 'toRotation');

end

