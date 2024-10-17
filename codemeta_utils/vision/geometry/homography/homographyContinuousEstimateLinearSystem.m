function [A,b]=homographyContinuousEstimateLinearSystem(x,dx)
xHom=homogeneous(x,3);
dxHom=homogeneous(dx,3,'velocities');

A=matStack(multikron(permute(xHom,[3 1 2]),hat3(xHom)));
b=matStack(multiprod(hat3(xHom),permute(dxHom,[1 3 2])));
