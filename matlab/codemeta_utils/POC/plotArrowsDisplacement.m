%function plotArrowsDisplacement(xStart,xDisplacement,varargin)
%Same as plot arrows, but give start and displacement instead of start and
%end of each arrow
function plotArrowsDisplacement(xStart,xDisplacement,varargin)
xEnd=xStart+xDisplacement;
plotArrows(xStart,xEnd,varargin{:})
