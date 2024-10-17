%function a=roterror(R,R1)
%Compute the maximum rotation error between two rotation matrices.
%The rotation error is quantified as the maximum angle between u=R*x and
%v=R1*x over all the possible x~=0
%
%EXAMPLE
%
%>> roterror(rot([0;0;1]*pi/6),eye(3))*180/pi
%ans =
%  30.0000
%
function a=roterror(R,R1)
a=acos(min(real(eig(R1'*R))));
