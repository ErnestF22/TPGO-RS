function [ c_dot ] = rotBundle_piDiff( c_bar, c_bar_dot )
% Maps the tangent vector c_bar_dot on TTSO(3) onto TSO(3) using the
% bundle map dpi(c_bar)
% See rotBundle_piDiffVec for the vector representation of the tangent vectors
% INPUTS:
%   c_bar := curve on the tangent bundle TSO(3) given as a [6 x 3] matrix,
%       [R;R*what]
%   c_bar_dot := tangent vector on TTSO(3) given as a [6 x 3] matrix,
%       [R*what; R*vhat]
% OUTPUTS:
%   c_dot := tangent vector of the curve c generating c_bar as a [3 x 3]
%      matrix, [R*what]

dpi = rotBundle_piMat(c_bar);

c_dot = dpi*c_bar_dot;
end

