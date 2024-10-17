%Returns the homography flow parameters from the motion and structure
%function alpha=homFlowParametersFromMotion(v,w,NVec)
%Inputs
%   v   [3 x 1] linear velocity
%   w   [3 x 1] angular velocity
%   n   [4 x 1] plane normal (in camera's frame)
%Output
%   alpha [8 x 1] vector with parameters of the infinitesimal homography
function alpha=homFlowParametersFromMotion(v,w,NVec)
n=planeNVecToNScaled(NVec);
alpha=[
    -w(2,:)+n(3,:)*v(1,:);      %1
    n(1,:)*v(1,:)-n(3,:)*v(3,:);%2
    n(2,:)*v(1,:)+w(3,:);       %3
    +w(1,:)+n(3,:)*v(2,:);      %4
    n(1,:)*v(2,:)-w(3,:);       %5
    n(2,:)*v(2,:)-n(3,:)*v(3,:);%6
    -n(2,:)*v(3,:)+w(1,:);      %7
    -w(2,:)-n(1,:)*v(3,:);      %8
    ];
