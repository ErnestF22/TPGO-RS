function [d,tMin]=essential_dist(Q1,Q2,varargin)
[tMin,d]=essential_distMinAngle(Q1,Q2,varargin{:});
d=sqrt(max(d,0));
