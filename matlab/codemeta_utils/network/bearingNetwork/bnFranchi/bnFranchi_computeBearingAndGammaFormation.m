function [y,g,E]=bnFranchi_computeBearingAndGammaFormation(x)
[d,NNodes]=size(x);
E=bnFranchi_graphEdges(NNodes);
y=bnFranchi_computeBearing(eye(d),x(:,E(:,1)),x(:,E(:,2)));
y12=y(:,1);
y21=y(:,1);
g=zeros(1,NNodes-2);
for iNode=1:NNodes-2
    y1i=y(:,2+2*iNode-1);
    y2i=y(:,2+2*iNode);
    g(iNode)=bnFranchi_Gamma(y2i,y12,y21,y1i);
end
