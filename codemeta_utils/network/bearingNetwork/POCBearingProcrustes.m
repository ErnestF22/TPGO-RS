function POCBearingProcrustes
global x0

flagDisplay=false;
resetRands
[xt,~,x0]=real_randGeodFun(rand(2,2));
N=50;
t=linspace(0,1,N);
xEval=funEvalVec(xt,t);
xAlign=funEvalVec(@align,xEval);
if flagDisplay
    figure(1)
    plotPairSequence(xEval,'b')
    hold on
    plotPairSequence(xAlign,'r')
    axis equal
end
figure(2)
plot(cnorm(reshape(xEval-xAlign,4,[])))

function xAlign=align(x)
global x0
xRef=cnormalize(x0);
[U,~,V]=svd(cnormalize(x)*xRef');
xAlign=U*V'*xRef;

function plotPairSequence(x,c)
N=size(x,3);
plotLines(zeros(2,N),squeeze(x(:,1,:)),c)
hold on
plotLines(zeros(2,N),squeeze(x(:,2,:)),c)
plotLines(squeeze(x(:,1,:)),squeeze(x(:,2,:)),c)
plotPoints(x(:,:,1),'rx')
plotPoints(x(:,:,end),'gx')
hold off
