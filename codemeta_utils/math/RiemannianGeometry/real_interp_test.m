function real_interp_test
NSamples=30;
t=linspace(0,3,NSamples);
tQuery=linspace(0,max(t),10*NSamples);
x=cos(2*pi*t)+t;
xQuery=real_interp(t,T,tQuery,'spline');
figure(1)

T=[cos(2*pi*t)+t; sin(2*pi*t+1); cos(2*pi*t+2)];
TQuery=real_interp(t,T,tQuery,'spline');
figure(2)
plot(t,T,'o',tQuery,TQuery,'-')

X0=randn(3,5);
R=repmat(eye(3),1,1,size(T,2));
X=rigidTransform(R,T,X0);
XQuery=real_interp(t,X,tQuery,'spline');
figure(3)
plotPoints(X,'bo')
hold on
for iPoint=1:size(XQuery,2)
    plotPoints(squeeze(XQuery(:,iPoint,:)),'b-')
end
hold off
