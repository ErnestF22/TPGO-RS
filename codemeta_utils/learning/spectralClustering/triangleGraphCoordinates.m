function xT=triangleGraphCoordinates(x,ET)
NTriangles=size(ET,1);
d=size(x,1);
xT=zeros(d,NTriangles);
for iTriangle=1:NTriangles
    xT(:,iTriangle)=mean(x(:,ET(iTriangle,:)),2);
end
