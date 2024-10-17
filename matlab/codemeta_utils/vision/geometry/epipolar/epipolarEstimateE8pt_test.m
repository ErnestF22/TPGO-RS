function epipolarEstimateE8pt_test
[G1,G2,x1,x2,X]=epipolar_dataset();
%testNetworkDisplay(cat(3,G1,G2),'Points',X,'references')
ETruth=epipolarBuildEFromG(G1,G2,'references');
EEst=epipolarEstimateE8pt(x1,x2);
disp([ETruth EEst])

