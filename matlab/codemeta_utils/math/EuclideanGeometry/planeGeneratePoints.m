function X=planeGeneratePoints(NVec,NPoints)
[N,d]=planeNVecToNd(NVec);
X=orthComplement(N)*randn(2,NPoints)+d*N*ones(1,NPoints);
