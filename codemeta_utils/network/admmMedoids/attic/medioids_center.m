function [mu,cMin]=medioids_center(x,varargin)
d=size(x,1);
mu=zeros(d,1);
if size(x,2)==0
    mu=NaN(d,1);
end

cMin=0;
for iDim=1:d
    c=medioids_center1dCost(x,'dim',iDim,varargin{:});
    [cMinDim,idxMin]=min(c);
    cMin=cMin+cMinDim;
    mu(iDim)=x(iDim,idxMin);
end
