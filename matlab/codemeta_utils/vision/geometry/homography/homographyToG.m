%Extract pose information and normal from a homography H
%function [G,N,lambda]=homographyToG(H,x1,x2)
%The image x1, x2 are used to solve the twisted pair ambiguity
function [G,N,lambda]=homographyToG(H,x1,x2)
[R,T,N,lambda]=homographyToRT(H,x1,x2);
G=RT2G(R,T);
