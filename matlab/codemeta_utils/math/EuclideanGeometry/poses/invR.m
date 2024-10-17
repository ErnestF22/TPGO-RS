%Compute the inverse (i.e., transpose) of a vector of rotations
%function R=invR(R)
function R=invR(R)
R=permute(R,[2 1 3]);
