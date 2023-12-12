function POCBearingDifferenceControl
X1=[0;0];
X2=[1;1];
x0=[2; 0];
%dx=[-1;0.1];
dx=cnormalize(randn(2,1));

Y1=@(x) cnormalize(x-X1);
Y2=@(x) cnormalize(x-X2);

Y1g=Y1(x0);
Y2g=Y2(x0);

c=@(x) Y1(x)'*Y2(x);
cg=Y1g'*Y2g;

xt=@(t) x0+t*dx;
disp(dx)

ct=@(t) c(xt(t));
t=linspace(0,4);
subplot(2,1,1)
plotfun(@(t) acos(ct(t)),t)
subplot(2,1,2)
xtt=evalfun(xt,t);
plot([X1(1) X2(1)],[X1(2) X2(2)],'r*')
hold on
plot(xtt(1,:),xtt(2,:))
hold off
axis([-1 3 -3 3])
axis square
axis equal
grid on
