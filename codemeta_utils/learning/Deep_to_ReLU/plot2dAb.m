% Plot a 2-D plane (line) from Ab_set
function plot2dAb(Ab_set, L, ax)
axXLim=ax(1,:);
axYLim=ax(2,:);
flagHold=ishold();
for i=1:L
    if i == 1
        Ai = Ab_set(1).A;
        bi = Ab_set(1).b;
    else
        Ai = Ab_set(i).A*Ab_set(i-1).A_cum;
        bi = Ab_set(i).A*Ab_set(i-1).b_cum+Ab_set(i-1).b;
    end
    [d,n] = size(Ai);
    for j=1:d
        vertices = zeros(2,0);
        p = Ai(j,:)';
        q = bi(j);
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
