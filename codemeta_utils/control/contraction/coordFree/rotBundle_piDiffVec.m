function [ c_dot_vhat ] = rotBundle_piDiffVec( c_bar, c_bar_dot_vhat )
% Maps the tangent vector c_bar_dot_vhat on TTSO(3) onto TSO(3) using the
% bundle map dpi(c_bar)
% See rotBundle_piDiff for the matrix representation of the tangent vectors
% INPUTS:
%   c_bar := curve on the tangent bundle TSO(3) given as a [6 x 3] matrix,
%       [R;R*what]
%   c_bar_dot_vhat := tangent vector on TTSO(3) given as a [6 x 1] vector,
%       [w;v]
% OUTPUTS:
%   c_dot_vhat := tangent vector of the curve c generating c_bar as a 
%    [3 x 1] vector

dpi = rotBundle_piDiffMat(c_bar);

c_dot_vhat = dpi*c_bar_dot_vhat;
end

