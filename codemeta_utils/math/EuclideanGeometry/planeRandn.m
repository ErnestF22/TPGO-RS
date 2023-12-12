%Generate a random plane
%function NVec=planeRandn(NVec,sigmaNormal,sigmaDistance)
function NVec=planeRandn(NVec,sigmaNormal,sigmaDistance)
if ~exist('sigmaNormal','var') || isempty(sigmaNormal)
    sigmaNormal=100;
end
if ~exist('sigmaDistance','var') || isempty(sigmaDistance)
    sigmaDistance=3;
end
    
[N,d]=planeNVecToNd(NVec);
N=sphere_randn(N,sigmaNormal);
d=d+sigmaDistance*randn(size(d));
NVec=planeNdToNVec(N,d);
