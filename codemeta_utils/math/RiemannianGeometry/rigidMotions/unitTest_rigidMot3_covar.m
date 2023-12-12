function [tests] = unitTest_rigidMot3_covar
% Test the covariant derivate on SE(3)
close all; clc;
tests = functiontests(localfunctions);
end

function testGeodesic(testCase)
% Test if 2 left-invar VF produces a valid tangent vector in se(3)
% [Gt,dGt,G0,dG0,vGVec,ddGt,dvGVec] = rigidMot3_randGeoFun();
xVec = randn(3,1);
yVec = randn(3,1);
xVel = randn(3,1);
yVel = randn(3,1);
X = @(R,d) [R*hat3(xVec) xVel;0 0 0 1];
Y = @(R,d) [R*hat3(yVec) yVel;0 0 0 1];

F1 = rigidMot3_covar(X,Y);
R = rot_randn; d = randn(3,1);
omega = rigidMot3_extractRot(F1(R,d));
v = rigidMot3_extractPos(F1(R,d));

R'*omega
v
end