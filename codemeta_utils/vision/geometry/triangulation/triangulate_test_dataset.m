function triangulate_test_dataset
%calibration='';
calibration='calibrated';

R=cat(3,eye(3),rot([0,1,0]*-pi/8));
T=cat(2,zeros(3,1), [3;1;0]);
G=RT2G(R,T);
if isempty(calibration)
    K=[320 0 320; 0 240 240; 0 0 1];
else
    K=eye(3);
end

X=[
     1  1 10;
     1 -1 10;
    -1  1 10;
    -1 -1 10;
     1  1 11;
     1 -1 11;
    -1  1 11;
    -1 -1 11]';

P=RTK2P(R,T,K);
x=projectFromP(P,X);
disp(x)

testNetworkDisplay(G,'Points',X)

save([mfilename '_data' calibration])
