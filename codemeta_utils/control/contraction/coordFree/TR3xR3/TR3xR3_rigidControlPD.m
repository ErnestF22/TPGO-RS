function [ dv ] = TR3xR3_rigidControlPD( x, vd, m, kd, kv, gradf_control )
%Controller for a rigid point mass in 3D with a reference trajectory. NOTE:
%the desired position has dynamics and new pos of x(7:9)
% INPUTS:
%   x := current state [9x1], [curr. pos.; curr. vel.;ref. pos]
%   vd := desired final velocity [3x1]
%   m := mass
%   kd, kv := control gains
%   gradf_control(x,xd) := function handle to gradient of cost function for
%       the controller

dv = m*(-kd*gradf_control(x(1:3),x(7:9)) - kv*(x(4:6)-vd));
end

