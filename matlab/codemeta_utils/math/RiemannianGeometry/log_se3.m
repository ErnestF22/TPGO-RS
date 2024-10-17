%Logarithm map SE(3)->se(3)
%
%function [w,v]=log_se3(e)
%
%Input
%   e   [4x4] transformation matrix
%Outputs
%   w   [3x1] vector with rotation axis times angle
%   v   [3x1] vector with translation part
%
%See also exp_se3

function [w,v]=log_se3(e)

if(e(4,1)~=0 && e(4,2)~=0 && e(4,3)~=0 && e(4,4)~=1)
    error('e is not in SE(3)');
end

R=e(1:3,1:3);
if(norm(R'*R-eye(3),'fro')>1e-6)
    error('e is not in SE(3), bad rotation');
end

p=e(1:3,4);

w=logrot(R);
theta=norm(w);

if(theta==0)
    v=p;
else
    wnorm=w/norm(w);
    
%     A=[R-eye(3);theta*wnorm'];
%     b=[hat(wnorm)*p; wnorm'*p];
%     
%     v=A\b;
    A=(eye(3)-R)*hat(wnorm)+wnorm*w';
%    v=A\p;
    v=A\(p*norm(w));
end
