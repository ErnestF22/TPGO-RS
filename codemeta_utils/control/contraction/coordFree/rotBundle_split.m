function [ p, u ] = rotBundle_split( Q )
% Split the point Q on TSO(3) into p in SO(3) and u in T_p SO(3)
% INPUTS:
%   Q := a point on TSO(3) represented by a [6 x 3] matrix
% OUTPUTS:
%   p := a point on SO(3) represented by a [3 x 3] matrix
%   q := a tangent vector on T_p SO(3) represented by a [3 x 3] matrix

p = Q(1:3,:);
u = Q(4:6,:);

end

