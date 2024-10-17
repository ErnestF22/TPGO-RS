function planeNScaledToNVec_test
NVec=randn(4,1);
n=planeNVecToNScaled(NVec);
NVecRec=planeNScaledToNVec(n);
subspace(NVec,NVecRec)