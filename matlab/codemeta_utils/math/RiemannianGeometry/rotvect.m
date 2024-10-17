%function R = rotvect(w)
%Same as rot() but creates an array of rotations using the columns of w
function R = rotvect(w)

temp=sqrt(sum(w.^2));
temp(temp<1e-12)=0;
theta=temp;
w=w./repmat(theta,3,1);
for(r=1:size(w,2))
    R(:,:,r) = eye(3,3) + hat(w(:,r))*sin(theta(r)) + hat(w(:,r))^2*(1-cos(theta(r)));
end
