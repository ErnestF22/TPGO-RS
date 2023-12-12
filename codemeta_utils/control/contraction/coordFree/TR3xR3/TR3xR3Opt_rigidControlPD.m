function [ dv ] = TR3xR3Opt_rigidControlPD( x, vd, m, gradf_control )
%Last Edited: Oct. 27 2020 by Bee Vang
% Controller for a rigid point mass in 3D with a reference trajectory. 
% In this controller, the gains are changing by solving an optimization problem.
% NOTE: the desired position has dynamics and new pos of x(7:9)
% INPUTS:
%   x := current state [12x1], [curr. pos.; curr. vel.;ref. pos;kd;kv;kref]
%   vd := desired final velocity [3x1]
%   m := mass
%   kd, kv := initial control gains (not used)
%   gradf_control(x,xd) := function handle to gradient of cost function for
%       the controller

% Extract the current gains
kd = x(10); kv = x(11);
% Compute the controller in the same way
dv = m*(-kd*gradf_control(x(1:3),x(7:9)) - kv*(x(4:6)-vd));
end

