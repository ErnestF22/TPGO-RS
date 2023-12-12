%function GPlane=planeToG(NVec,varargin)
%Return a reference frame with the origin in the plane and the z axis
%aligned to the normal of the plane.
function GPlane=planeToG(NVec,varargin)
[RPlane,TPlane]=planeToRT(NVec,varargin{:});
GPlane=RT2G(RPlane,TPlane);
