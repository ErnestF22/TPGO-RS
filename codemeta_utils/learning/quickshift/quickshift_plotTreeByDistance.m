function quickshift_plotTreeByDistance(treeVectorMembership,treeVectorDistance,treeVectorDensity)
NPoints=length(treeVectorMembership);

if exist('treeVectorDensity','var')
    y=treeVectorDensity;
else
    y=1:NPoints;
end

%vector of final locations
x=zeros(1,NPoints);

xPrev=ones(1,NPoints);
while ~isequal(x,xPrev)
    xPrev=x;
    x=x(treeVectorMembership)+treeVectorDistance;
end

plotLines([x;zeros(1,NPoints)],[x;y],'b')
hold on
plotLines([x;y],[x(treeVectorMembership);y],'b')
hold off
for iPoint=1:NPoints
    text(x(iPoint),0,num2str(iPoint));
end
