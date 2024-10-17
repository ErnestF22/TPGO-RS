%Plot field given by position and direction vectors
%function plotField(x,dx,varargin)
%   x   [d x N] array with positions
%   dx  [d x N] array with directions
function plotField(x,dx,varargin)
flagHold=ishold();

plotPoints(x)
hold on

d=size(x,1);
switch d
    case 2
        quiver(squeeze(x(1,:,:)),squeeze(x(2,:,:)),squeeze(dx(1,:,:)),squeeze(dx(2,:,:)),varargin{:})
    case 3
        quiver3(squeeze(x(1,:,:)),squeeze(x(2,:,:)),squeeze(x(3,:,:)),squeeze(dx(1,:,:)),squeeze(dx(2,:,:)),squeeze(dx(3,:,:)),varargin{:})
    otherwise
        error('First dimension of the data must be two or three')
end
        
if ~flagHold
    hold off
end
