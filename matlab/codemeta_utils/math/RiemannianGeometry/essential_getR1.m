%Get first rotation of a point in the QREM
%function R=essential_getR1(Q)
%Input
%   Q   [6x3] matrix, a point in the QREM
%Output
%   R   [3x3] rotation, the first rotation in Q
function R=essential_getR1(Q)
R=Q(1:3,:,:);
