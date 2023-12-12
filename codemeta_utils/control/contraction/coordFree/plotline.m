%function plotline(x,y,n,c,style)
%Plots a line in the 2D region [min(x), max(x)] X [min(y), max(y)]
%satisfying n'*[x;y]+c=0 
function plotline(x,y,n,c,style)
if ~exist('style','var')
    style='-';
end


%coordinates for the four corners
minx=min(x);
miny=min(y);
maxx=max(x);
maxy=max(y);

%sign of line function
f=@(x,y) sign(n'*[x;y]+c);

%evaluate sign of line function at the four corners of the region
fll=f(minx,miny);
flh=f(minx,maxy);
fhh=f(maxx,maxy);
fhl=f(maxx,miny);

%make sure line is visible
if ~(fll==flh && flh==fhh && fhh==fhl)
    %find to points of the line as the intersection with the boundaries

    p=1; %counts which intersection we are looking for (first or second)
    if fll~=fhl     %south
        yp(p)=miny;
        %solve n(1)x+n(2)y+c=0
        xp(p)=(-n(2)*yp(p)-c)/n(1);
        p=p+1;
    end
    if fll~=flh     %west
        xp(p)=minx;
        yp(p)=(-n(1)*xp(p)-c)/n(2);
        p=p+1;
    end
    if fhl~=fhh     %east
        xp(p)=maxx;
        yp(p)=(-n(1)*xp(p)-c)/n(2);
        p=p+1;
    end
    if flh~=fhh     %north
        yp(p)=maxy;
        xp(p)=(-n(2)*yp(p)-c)/n(1);
        p=p+1;
    end

    plot(xp,yp,style,'MarkerSize',3,'LineWidth',3)
    axis([minx maxx miny maxy])
end

