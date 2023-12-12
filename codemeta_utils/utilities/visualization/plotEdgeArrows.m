%Plot lines between points given by edge list
%function plotEdges(x,E,varargin)
%Inputs
%   x   [d x NPoints] coordinates of the points
%   E   [NEdges x 2] edge list
%Optional arguments are passed directly to plotLines
function plotEdgeArrows(x,E,varargin)
xStart=x(:,E(:,1));
xEnd=x(:,E(:,2));
plotArrows(xStart,xEnd,varargin{:})
