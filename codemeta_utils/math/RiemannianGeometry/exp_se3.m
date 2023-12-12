%Exponential map se(3)->SE(3)
%
%function e=exp_se3(w,v)
%
%Inputs
%   w   [3x1] vector with rotation axis times angle
%   v   [3x1] vector with translation part
%Output
%   e   [4x4] transformation matrix
%
%See also exp_se3

%%AUTORIGHTS%%

function e=exp_se3(w,v)
theta=norm(w);

if(theta>1e-15)
%     S=[hat(w)/theta v; zeros(1,4)];
%     e=eye(4)+theta*S+(1-cos(theta))*S^2+(theta-sin(theta))*S^3;
   wnorm=w/theta;

    R=rot(w);
    p=(eye(3)-R)*hat(wnorm)*v+w'*v*wnorm;
    
    e=[R p/norm(w); zeros(1,3) 1];
else
    e=[eye(3) v; zeros(1,3) 1];
end