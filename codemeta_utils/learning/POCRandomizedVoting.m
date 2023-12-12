function POCRandomizedVoting
resetRands(2)

NPlanes=4;
b=squeeze(sphere_randn(eye(2,1),[],NPlanes));
d=randn(1,NPlanes);
bHom=[b;d];

for iPlane=1:NPlanes
    draw2dLine(bHom(:,iPlane))
    hold on
end
hold off

NPointsGroup=round(20*(rand(1,NPlanes)+0.5));

XIn=[];
for iPlane=1:NPlanes
    preX=1.5*randn(2,NPointsGroup(iPlane));
    XIn=[XIn homogeneousProjectPoints(preX,bHom(:,iPlane))];
end

NPointsIn=size(XIn,2);
XOut=1.5*randn(2,round(NPointsIn/5));
X=[XIn XOut];

figure(1)
hold on
plotPoints(XIn)
plotPoints(XOut,'r')
hold off

NPoints=size(X,2);
A=zeros(NPoints,NPoints);
threshold=0.01;

flagDebug=false;

NSamples=10000;
for iSample=1:NSamples
    idxSample=randperm(NPoints,2);
    bOrth=fitLine(X(:,idxSample));
    r=residuals(X,bOrth);
    idxInliers=r<threshold;
    if flagDebug
        figure(4)
        plotPoints(X,'b')
        hold on
        plotPoints(X(:,idxInliers),'g')
        plotPoints(X(:,idxSample),'r')
        hold off
        pause
    end
        
    A(idxInliers,idxInliers)=A(idxInliers,idxInliers)+1;
end

figure(2)
imagesc(A)

figure(3)
plot(sum(A))

function r=residuals(x,bOrth)
r=sum((x-homogeneousProjectPoints(x,bOrth)).^2);

function bOrth=fitLine(x)
[U,~,~]=svd(homogeneous(x,3));
bOrth=U(:,3);

