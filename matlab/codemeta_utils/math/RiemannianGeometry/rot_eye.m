%function R=rot_eye(R)
%Returns the identity for the space of rotations.
%If R is a matrix, returns eye(size(R))
%If R is a scalar, returns eye(R)
%If R is omitted, returns eye(3)
function R=rot_eye(R)
if ~exist('R','var')
    R=eye(3);
else
    if isscalar(R)
        R=eye(R);
    else
        R=repmat(eye(size(R,1)),1,1,size(R,3));
    end
end
