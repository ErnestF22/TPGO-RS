%Get second rotation of a point in the QREM
%function R=essential_getR2(Q)
%Input
%   Q   [6x3] matrix, a point in the QREM
%Output
%   R   [3x3] rotation, the second rotation in Q
function R=essential_getR2(Q)
R=Q(4:6,:,:);
