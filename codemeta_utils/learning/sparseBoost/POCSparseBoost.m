function POCSparseBoost
resetRands()
d=2;
NPointsPerClass=50;
NLearners=1000;

L=3;

switch 1
    case 1
        X1=randn(d,NPointsPerClass);
        y1=ones(1,NPointsPerClass);
        X2a=randn(d,NPointsPerClass)+L;
        X2b=X2a;
        X2c=X2a;
        X2b(1,:)=X2b(1,:)-L;
        X2c(2,:)=X2c(2,:)-L;
        X2=[X2a X2b X2c];
        y2=-ones(1,3*NPointsPerClass);
    case 2
        XPre=1*randn(d,4*NPointsPerClass);
        nXPre=sqrt(sum(XPre.^2));
        fClass=nXPre>1;
        X1=XPre(:,fClass);
        y1=ones(1,sum(fClass));
        X2=XPre(:,~fClass);
        y2=-ones(1,sum(~fClass));
end
X=[X1 X2];
y=[y1 y2];
NPoints=size(X,2);

switch 2
    case 1
        NLearnersHalf=round(NLearners/2);
        wLearners=[ones(1,NLearnersHalf); zeros(1,NLearnersHalf)];
        wLearners=[wLearners flipud(wLearners)];
        bLearners=linspace(-4,4,NLearnersHalf);
        bLearners=[bLearners bLearners];
    case 2
        wLearners=cnormalize(randn(d,NLearners));
        bLearners=10*(rand(1,NLearners)-0.5);
end

rLearners=max(0,wLearners'*X+bLearners'*ones(1,NPoints));
rLearners=[rLearners; max(0,-(wLearners'*X+bLearners'*ones(1,NPoints)))];

C=1;

u=ones(1,NPoints);
cvx_begin
    variables w(2*NLearners,1) b z(NPoints,1)
    minimize (norm(w,1)+C*sum(z))
    subject to
        y.*(w'*rLearners+b)-u+z'>=0
        z>=0
cvx_end


%w(100)=1;
xGrid=linspace(-3,4);
figure(1)
z=imagefun(xGrid,xGrid,@(x) fEnsamble(flipud(x),w,b,wLearners,bLearners));
colormap gray
hold on
plotPoints(X1)
plotPoints(X2,'r')
hold off

figure(2)
contour(xGrid,xGrid,z,[0 0])
hold on
plotPoints(X1)
plotPoints(X2,'r')
hold off

figure(3)
plot(sign(fEnsamble(X,w,b,wLearners,bLearners)))

disp(sum(abs(w)>1e-2))


function f=fSVM(x,w,b)
f=w'*x+b;

function f=fEnsamble(x,w,b,wLearners,bLearners)
NPoints=size(x,2);
rLearners=max(0,wLearners'*x+bLearners'*ones(1,NPoints));
rLearners=[rLearners; max(0,-(wLearners'*x+bLearners'*ones(1,NPoints)))];

f=w'*rLearners+b;
