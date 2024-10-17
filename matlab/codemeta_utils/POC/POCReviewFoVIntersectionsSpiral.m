function POCReviewFoVIntersectionsSpiral

NPoints=30;
t=linspace(0,2*pi,NPoints);
V=3.5;
sigma=V/2;
Rout=15;
Rin=5;
cs=spline([t(1) t(end)],[0.5 t(1) t(end) 1.5]);
theta=ppval(cs,t);
figure(2)
plot(t,theta)

r=Rout-(Rout-Rin)*(t-min(t))/(max(t)-min(t));

p=multiprod([cos(theta);sin(theta)],r,1,1);

figure(1)
plotPoints(p)
hold on

for iPoint=1:NPoints
    pThis=p(:,iPoint);
    d=sqrt(euclideanDistMatrix(pThis,p));
    flagNeighbors=d<V;
    flagNeighbors(iPoint)=false;
    NNeighbors=sum(flagNeighbors);
    pNeighbors=p(:,flagNeighbors);
    plotLines(repmat(pThis,1,NNeighbors),pNeighbors);
    for iNeighbor=1:NNeighbors
        U=cnormalize(pNeighbors(:,iNeighbor)-pThis);
        c=pThis+sigma*U;
        %plotPoints(c,'x')
        plotCircle(c(1),c(2),sigma);
    end
end

hold off
axis equal
