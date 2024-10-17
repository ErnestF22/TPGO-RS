function POCSVM
resetRands()
d=2;
NPointsPerClass=30;
X1=randn(d,NPointsPerClass);
y1=ones(1,NPointsPerClass);
X2=randn(d,NPointsPerClass)+2;
y2=-ones(1,NPointsPerClass);
X=[X1 X2];
y=[y1 y2];
NPoints=size(X,2);


C=1;

u=ones(1,NPoints);
cvx_begin
    variables w(d,1) b z(NPoints,1)
    minimize (norm(w,1)+C*sum(z))
    subject to
        y.*(w'*X+b)-u+z'>=0
        z>=0
cvx_end

xGrid=linspace(-3,4);
figure(1)
z=imagefun(xGrid,xGrid,@(x) fSVM(x,w,b));
colormap gray
hold on
plotPoints(X1)
plotPoints(X2,{'Color','r'})
hold off

figure(2)
plot(sign(fSVM(X,w,b)))
function f=fSVM(x,w,b)
f=w'*x+b;
