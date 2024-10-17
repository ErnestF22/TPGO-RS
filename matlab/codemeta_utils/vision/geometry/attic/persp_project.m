%PERSP_PROJECT		x = persp_project(q,f)
%	Computes the perspective projection x=f*q/q(3) of a set of points.
%	Here both q and x are 3 by n by m vectors where n is the number 
%       of points, m is the number of frames and f is the focal length 
%       of the camera.
%
%	See also sphere_project

function x = persp_project(q,f)

if nargin == 1
   f=1;
end

%x = f*q./(ones(3,1)*q(3,:));
x = f*q./repmat(q(3,:,:),3,1);
