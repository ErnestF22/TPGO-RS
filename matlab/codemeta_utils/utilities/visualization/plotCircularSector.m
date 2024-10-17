%Plot a 2-D circular sector with a given center
%function plotCircularSector(o,tStart,tEnd,radius,color,varargin)
%Inputs
%   o               coordinates of the center
%   tStart, tEnd    starting and ending absolute angles
%   radius          radius
%   color           filling color
%All optional parameters are passed straight to the patch function
function plotCircularSector(o,tStart,tEnd,radius,varargin)
if ~exist('radius','var') || isempty(radius)
    radius=1;
end
if isempty(varargin)
    varargin={'r'};
end

t=linspace(tStart,tEnd,round(abs(tEnd-tStart)/0.04));
x=[zeros(2,1) radius*[cos(t); sin(t)]];

if ~ishold()
    cla;
end
patch(x(1,:)+o(1),x(2,:)+o(2),varargin{1},'edgeColor','none',varargin{2:end})
