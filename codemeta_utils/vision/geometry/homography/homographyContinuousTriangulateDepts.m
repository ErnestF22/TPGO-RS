function [lambda,xHom]=homographyContinuousTriangulateDepts(x,nVec)
xHom=homogeneous(x,3);
lambda=-1./(nVec'*xHom);
