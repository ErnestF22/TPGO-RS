%function [N,d]=planeNVecToNd(NVec)
%Pass from the plane representation NVec'*[X;1]=0 to the representation
%N'*x=d
function [N,d]=planeNVecToNd(NVec)
NVec=planeNormalize(NVec);
N=NVec(1:3,:);
d=-NVec(4,:);