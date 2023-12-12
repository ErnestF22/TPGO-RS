function [ dpi ] = rotBundle_piDiffMat( c_bar )
% Return the differential of the bundle map, pi, as a matrix evaluated at 
% c_bar (a point on TSO(3)) (For our case, dpi is a const [3 x 6] matrix)
% NOTE: due to our choice of representation for a point on the tangent
% bundle TSO(3), piMat (pi) and piDiffMat (dpi) are the same

dpi = [eye(3), zeros(3,3)];

end

