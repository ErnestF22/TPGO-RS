%function y=sphere_eye(y)
%Returns the identity for the sphere.
%If y is a matrix, returns eye(size(y))
%If y is a scalar, returns eye(y)
%If y is omitted, returns eye(3,1)
function y=sphere_eye(y)
if ~exist('y','var')
    y=eye(3,1);
else
    if isscalar(y)
        y=eye(y,1);
    else
        y=eye(size(y));
    end
end
