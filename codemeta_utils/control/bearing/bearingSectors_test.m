function bearingSectors_test
resetRands()
%%
NNodes=10;
x=[zeros(2,1) 2*randn(2,NNodes-1)];
E=adj2edges(ones(NNodes));
r=1;

[y,a]=bearingNetworkCollisionImages(x,E,r);
[y1,a1]=bearingNetworkCollisionCollectImagesNode(y,a,E,1);
x1=x(:,1);
ay1=atan2(y1(2,:),y1(1,:));
%%

ag=deg2rad(-40);
g=5*[cos(ag); sin(ag)];
plotPoints(x)
hold on
bearingPlot(x1,y1)
for jNode=1:size(y1,2)
    plotCircularSector(x1,ay1(jNode)-a1(jNode),ay1(jNode)+a1(jNode),r);
end
bearingPlot(x1,g,'b')
[uProj,flagHit,dProj]=bearingSectorsFind(g,y1,a1);
display(flagHit)
display(dProj)
if isempty(uProj)
    disp('No projection!')
else
    bearingPlot(x1,uProj,'g')
    gProj=bearingSectorsProject(g,y1,a1);
    bearingPlot(x1,gProj,'c')
end
hold off
axis equal
%%

function y=angle2Vector(a)
y=[cos(a); sin(a)];

function a=vector2Angle(y)
a=atan2(y(2,:),y(1,:));

%%
function [yMerged,aMerged]=mergeSectors2D(y,a)
NNodes=size(y,2);
ay=atan2(y(2,:),y(1,:));

%compute endpoints of intervals
ayMax=modAngle(ay+a);
ayMin=modAngle(ay-a);

%sort data by min angle
[ayMin,idxSort]=sort(ayMin);
ayMax=ayMax(idxSort);
y=y(:,idxSort);
a=a(idxSort);

%do the merging
yMerged=y(:,1);
aMerged=a(1);
ayMergedMin=ayMin(1);
ayMergedMax=ayMax(1);

for iNode=2:NNodes
    ayiMax=ay(iNode)+a(iNode);
    ayiMin=ay(iNode)-a(iNode);
    
    if angleDiff(ayiMin,ayMergedMax(end))>0
    end
end


%%