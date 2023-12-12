%Plot a circular sector representing a field of view
%function bearingPlotFov(x,R,fov,color,y0)
%Inputs
%   x       camera center
%   R       2-D camera rotation (can also be expressed as absolute heading
%           angles or heading directions)
%   fov     field of view in rads. If it has two elements, the smaller is
%           used to overlap a second, smaller, white circular sector
%   color   color to use to fill in the sector (default, 'r')
%   y0      direction of the center of the f.o.v. in local coordinates
%           (default: [1;0])
function bearingPlotFov(x,R,fov,color,y0)
if ~exist('y0','var') || isempty(y0)
    y0=[1;0];
end
if ~exist('color','var') || isempty(color)
    color='r';
end
radius=1;

N=size(x,2);

if size(R,1)==1 || (size(R,1)==2 && size(R,2)==N && size(R,3)==1)
    R=bearingHeading2Rotation(R);
end

if N>1
    flagHold=ishold();
    for iN=1:N
        bearingPlotFov(x(:,iN),R(:,:,iN),fov,color,y0)
        hold on
    end
    if ~flagHold
        hold off
    end
else
    Ry0=R*y0;
    t0=atan2(Ry0(2),Ry0(1));
    fovMax=max(fov);
    plotCircularSector(x,t0-fovMax/2,t0+fovMax/2,radius,color);
    if length(fov)>1
       fovMin=min(fov); 
        plotCircularSector(x,t0-fovMin/2,t0+fovMin/2,radius*0.7,'w');
    end
end
