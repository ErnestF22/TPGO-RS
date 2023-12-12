function POCGeodesicSphereColinear
%study of conditions under which geodesics in the shape space S_2^3 pass
%through degenerate configurations (collinear vectors)

switch 2
    case 1
        resetRands(2)

        X0=sphere_randn(eye(4,1));
        Ra=householderRotation(X0(1:2),1);
        X=blkdiag(Ra,Ra)*X0;
        YColin=randSphereColinear();
        logXYColin=sphere_log(X,YColin);
        distXYColin=norm(logXYColin);
        distXY=(distXYColin+rand*(pi-distXYColin));
        logXY=distXY*cnormalize(logXYColin);
        Y=sphere_exp(X,logXY);
    case 2
        X=cnormalize([0.3; 0; 0.6; -0.11]);
        Y=cnormalize([-0.42; -0.25; -0.86; 0]);
end

logXY=sphere_log(X,Y);
Rx=diag([1;-1]);
RxY=blkdiag(Rx,Rx)*Y;
logXRxY=sphere_log(X,RxY);

gammaXY=@(t) sphere_exp(X,t*logXY);
gammaXY1=@(t) [eye(2) zeros(2)]*gammaXY(t);
gammaXY2=@(t) [zeros(2) eye(2)]*gammaXY(t);
gammaXRxY=@(t) sphere_exp(X,t*logXRxY);
gammaXRxY1=@(t) [eye(2) zeros(2)]*gammaXRxY(t);
gammaXRxY2=@(t) [zeros(2) eye(2)]*gammaXRxY(t);

% funPlot(gamma1)
% hold on
% funPlot(gamma2,[],'--')
% hold off

t=linspace(0,1,500);
figure(1)
subplot(2,1,1)
funPlot(@(t) subspace(gammaXY1(t),gammaXY2(t)),t)
subplot(2,1,2)
funPlot(@(t) subspace(gammaXRxY1(t),gammaXRxY2(t)),t)


gammaXY1Eval=funEval(gammaXY1,t);
gammaXY2Eval=funEval(gammaXY2,t);
gammaXRxY1Eval=funEval(gammaXRxY1,t);
gammaXRxY2Eval=funEval(gammaXRxY2,t);

figure(2)
quiver([0;0],[0;0],Y([1 3]), Y([2 4]),1,'r')
hold on
plot(gammaXY1Eval(1,:),gammaXY1Eval(2,:),'r--')
plot(gammaXY2Eval(1,:),gammaXY2Eval(2,:),'r--')
quiver([0;0],[0;0],RxY([1 3]), RxY([2 4]),1,'b')
plot(gammaXRxY1Eval(1,:),gammaXRxY1Eval(2,:),'b--')
plot(gammaXRxY2Eval(1,:),gammaXRxY2Eval(2,:),'b--')
quiver([0;0],[0;0],X([1 3]), X([2 4]),1,'g')

hold off
axis equal

t01=0.505;
subspace(gammaXY1(t01),gammaXY2(t01))
t02=0.524;
subspace(gammaXRxY1(t02),gammaXRxY2(t02))

keyboard

function Y=randSphereColinear()
k=sphere_randn(eye(2,1));
Y=kron(k,sphere_randn(eye(2,1)));
