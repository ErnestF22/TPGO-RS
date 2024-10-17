%Interpolate a 2-D surface from scattered samples on a square grid
%function [X,Y,Z]=interpGridScattered(x,y,z,NPoints)
%Given a 2-D function with values z(i) at points (x(i),y(i)), interpolate
%these samples onto a square grid of NPoint^2 points.
function [X,Y,Z]=interpGridScattered(x,y,z,NPoints)
xlin = linspace(min(x),max(x),NPoints);
ylin = linspace(min(y),max(y),NPoints);
[X,Y] = meshgrid(xlin,ylin);
f = scatteredInterpolant(shiftdim(x),shiftdim(y),shiftdim(z));
Z = f(X,Y);
