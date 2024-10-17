%Pass from homogeneous to scaled vector representation of a plane
%function NScaled=planeNVecToNScaled(NVec)
%The normal is scaled by the inverse of the distance.
%Note: this representation does not work for planes passing through the origin.
function NScaled=planeNVecToNScaled(NVec)
NScaled=planeNormalize(NVec,'unitDistance');
NScaled=NScaled(1:end-1,:);

