%Creates two views and points visible in the two views
%function [G1,G2,x1,x2,X]=epipolar_dataset()
%G1 and G2 are in the 'reference' interpretation
function [G1,G2,x1,x2,X]=epipolar_dataset()
[G,X]=testNetworkCreateAbsolutePoses(7);
G=invg(G(:,:,1:2));
G1=G(:,:,1);
G2=G(:,:,2);
x=projectFromG(G,X,'references');
x1=x(:,:,1);
x2=x(:,:,2);
