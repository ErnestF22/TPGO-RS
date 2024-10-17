function projectFromG_test
[G,X]=testNetworkCreateAbsolutePoses(7);
G1=G(:,:,1);
x=projectFromG(G1,X,'references');

xTruth=[eye(2) zeros(2,1)]*homogeneous([eye(3) zeros(3,1)]*invg(G1)*homogeneous(X,4),3);

disp(x-xTruth)
