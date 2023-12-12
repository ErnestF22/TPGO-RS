function [Gamma,thrust]=quadrotor_inputUnpack(dx)
Gamma=dx(1:3,:);
thrust=dx(4,:);

