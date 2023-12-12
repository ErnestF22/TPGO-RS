function bearingSectorsPlot(x,y,a,r)
if ~exist('r','var')
    r=0.5;
end
ay=atan2(y(2,:),y(1,:));
for jNode=1:size(y,2)
    plotCircularSector(x,ay(jNode)-a(jNode),ay(jNode)+a(jNode),r);
end
axis equal

