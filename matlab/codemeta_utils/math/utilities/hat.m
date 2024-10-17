%h = hat(p) computes the wedge operator for the cross product, that is
%hat(p)*q = cross(p,q). p is a 3 by n vector. h is a 3 by 3 by n matrix.
%
%See also vee

function h = hat(p)
warning('This function is deprecated, use hat3 instead')
n = size(p,2);
p=shiftdim(p,-1);
z = zeros(1,1,n);
h=[ z         p(1,3,:) -p(1,2,:);
   -p(1,3,:)  z         p(1,1,:);
    p(1,2,:) -p(1,1,:)  z];
