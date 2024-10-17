%Plot lines between points given by edge list
%function plotEdges(x,E,varargin)
%Inputs
%   x   [d x NPoints] coordinates of the points
%   E   [NEdges x 2] edge list
%Optional arguments are passed directly to plotLines
function plotEdges(x,E,varargin)
xStart=x(:,E(:,1));
xEnd=x(:,E(:,2));
plotLines(xStart,xEnd,varargin{:})
