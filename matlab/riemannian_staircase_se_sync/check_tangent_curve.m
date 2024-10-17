function is_tangent = check_tangent_curve(R0, v0)
%CHECK_ Summary of this function goes here
%   Detailed explanation goes here


% v0=rot_randTangentNormVector(R0);
% ct=@(t) retraction_stiefel(R0,t*v0); %retractions
h=1e-9; 
v0Approx = (retraction_stiefel(R0,h*v0)-R0)/h;

% v0Approx=funApproxDer(ct,0,1e-9);
disp([v0 v0Approx])


if max(abs(v0 - v0Approx), [], "all") < 1e-5
    is_tangent = boolean(1);
else
    is_tangent = boolean(0);

end

