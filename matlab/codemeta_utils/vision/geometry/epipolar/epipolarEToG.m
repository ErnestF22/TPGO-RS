function [G,lambda]=epipolarEToG(E,x1,x2)
[R,T,lambda]=epipolarEToRT(E,x1,x2);
G=RT2G(R,T);
