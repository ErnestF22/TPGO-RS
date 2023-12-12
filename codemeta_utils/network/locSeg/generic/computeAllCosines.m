function [c,v]=computeAllCosines(v)
v=cnormalize(v);
c=v'*v;
