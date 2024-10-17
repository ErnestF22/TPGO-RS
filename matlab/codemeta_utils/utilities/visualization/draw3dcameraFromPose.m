%function draw3dcameraFromPose(R,T,varargin)
%Draw a pyramid and a set of axes to represent a camera
%R  axes of the camera
%T  center of the camera
%Optional arguments: same as draw3dcameraFromAxesAndCenter
%
%See also draw3dcameraFromAxesAndCenter

%%AUTORIGHTS%%

function draw3dcameraFromPose(R,T,varargin)
if(nargin<1)
    R=eye(3);
end
if(nargin<2)
    T=[0;0;0];
end

%get axes and camera center from pose
T=-R'*T;
R=R';

draw3dcameraFromAxesAndCenter(R,T,varargin{:});
