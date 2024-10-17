%function H=homographyFromG(G,NVec)
%Compute homography from pose and plane. The plane vector NVec is a [4 x 1]
%vector and points belonging to the plane satisfy NVec'*[X;1]=0
function H=homographyFromG(G,NVec)
[R,T]=G2RT(G);
H=homographyFromRT(R,T,NVec);
