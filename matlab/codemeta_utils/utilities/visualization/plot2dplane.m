% Plot a 2-D plane (line) from its equation clipping againts axes
% function plot2dplane(n,d)
% Plot a line satisfying p'*x+q=0 against current axes
function plot2dplane(p,q)
vertices=zeros(2,0);
ax=axis;
axXLim=ax(1:2);
axYLim=ax(3:4);
flagHold=ishold();
for iLim=1:2
    y=intersectVertical(p,q,axXLim(iLim));
    if isInInterval(y,axYLim)
        vertices=[vertices [axXLim(iLim);y]];
    end
    x=intersectHorizontal(p,q,axYLim(iLim));
    if isInInterval(x,axXLim)
        vertices=[vertices [x;axYLim(iLim)]];
    end
end

if size(vertices,2)>0
    hold on
    %remove extra vertices when the plot goes through a corner of the axes
    vertices=uniquetol(vertices',1e-6,'ByRows',true)';
    plotPoints(vertices,'-')
    if ~flagHold
        hold off
    end
end

function y=intersectVertical(p,q,x)
%Intersect with the vertical line having x=constant
%Returns Inf if the lines are parallel
if p(2)==0
    y=Inf;
else
    y=-(p(1)'*x+q)/p(2);
end

function x=intersectHorizontal(p,q,y)
%Intersect with the vertical line having x=constant
%Returns Inf if the lines are parallel

%Implemented by transposing x/y axes
x=intersectVertical(flipud(p),q,y);

function flag=isInInterval(x,interval)
%Return true if x is between min(interval) and max(interval)
flag=and(x>=min(interval), x<=max(interval));
