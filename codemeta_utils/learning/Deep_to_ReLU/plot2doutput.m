% Plot a 2-D plane (line) given A and b -> Ax+b = 0
function plot2doutput(A_out,b_out,ax)
axXLim=ax(1,:);
axYLim=ax(2,:);
flagHold=ishold();
d = size(A_out,1);
for j=1:d
    vertices = zeros(2,0);
    p = A_out(j,:)';
    q = b_out(j);
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
