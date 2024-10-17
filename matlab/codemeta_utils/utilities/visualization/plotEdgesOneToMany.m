%Plot edges from a single origin to multiple endpoints
function plotEdgesOneToMany(xStart,xEnd,varargin)
Nx=size(xEnd,2);
E=[ones(1,Nx);(1:Nx)+1]';
plotEdges([xStart xEnd],E,varargin{:})
