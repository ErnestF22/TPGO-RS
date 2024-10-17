function epipolarEToRT_test
[G1,G2,x1,x2,X]=epipolar_dataset();
%testNetworkDisplay(cat(3,G1,G2),'Points',X,'references')
ETruth=epipolarBuildEFromG(G1,G2,'references');
G12=computeRelativePoseFromG(G1,G2,'references');
G12(1:3,4)=cnormalize(G12(1:3,4));
[R12,T12,lambda]=epipolarEToRT(ETruth,x1,x2);
disp([G12 RT2G(R12,T12)])
