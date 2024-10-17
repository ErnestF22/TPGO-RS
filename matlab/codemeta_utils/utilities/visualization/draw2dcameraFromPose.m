%function draw2dcameraFromPose(R,T,varargin)
%Draw a pyramid and a set of axes to represent a camera
%   R,T     pose of the camera (world to camera transformation)
%Optional arguments: same as draw3dcameraFromAxesAndCenter
%
%See also draw2dcameraFromAxesAndCenter

%%AUTORIGHTS%%

function draw2dcameraFromPose(R,T,varargin)
if(nargin<1)
    R=eye(3);
end
if(nargin<2)
    T=[0;0;0];
end

%get axes and camera center from pose
[R,T]=invRT(R,T);

draw2dcameraFromAxesAndCenter(R,T,varargin{:});
