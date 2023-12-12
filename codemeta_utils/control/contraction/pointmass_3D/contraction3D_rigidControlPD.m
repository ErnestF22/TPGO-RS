function [ dv ] = contraction3D_rigidControlPD( x, xd, vd, m, kd, kv, gradf_control )
%Controller for a rigid point mass in 3D
if exist('gradf_control','var')
    dv = m*(-kd*gradf_control(x(1:3),xd) - kv*(x(4:6)-vd));
else
    dv = m*(kd*(xd-x(1:3))/sqrt(1+norm(xd-x(1:3))^2/2) + kv*(vd-x(4:end)));
    % dv = m*(kd*(xd-x(1:3)) + kv*(vd-x(4:end))); %using p(x) = norm(x)^2
end
end

