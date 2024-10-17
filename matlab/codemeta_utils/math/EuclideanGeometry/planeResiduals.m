%Compute residuals of the plane equation
%function e=planeResiduals(NVec,X)
%The function returns the result of e=N'*X-d. For points X belonging to the
%plane, e should be zero
function e=planeResiduals(NVec,X)
[N,d]=planeNVecToNd(NVec);
e=N'*X-d;

