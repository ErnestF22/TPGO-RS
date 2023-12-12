function [ dx ] = TR3xR3_rigidModel( x, xd, m, F, gradf_ref, kp )
%Create a rigid point mass model in R^3 with a reference trajectoryy. 
%Model is given by m*ddx=F
% INPUT:
%   x := current state given as [9x1] vector
%   xd := desired final pos [3x1]
%   m := mass
%   F := control input
%   gradf_ref(x,xd) := function handle to gradient of cost function for
%       reference trajectory
%   kp := positive scalar gain

dx = [zeros(3) eye(3) zeros(3);zeros(3,9);zeros(3,9)]*x ... % state
    + [zeros(3);eye(3);zeros(3)]*F*1/m ... % control
    + [zeros(6,1);-kp*gradf_ref(x(7:9),xd)]; % reference traj.
end

