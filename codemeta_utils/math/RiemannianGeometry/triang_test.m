function triang_test
X=[1;1;-1.5];
x=triang_projectX(X);
XRec=triang_getX(x);
T=triang_tangentBasis(x);
f=@(z) triang_getX(triang_exp(x,reshape(T*z,2,2)));


figure(1)
draw3dcameraFromRT(eye(3),zeros(3,1))
hold on
draw3dcameraFromRT(eye(3),[0;0;-1],'references')
plotPoints(X)
plotPoints([x; 1 0],{'Color',[1 0 0]})
plotSphereMap(f,'stepRadii',0.1,'NRadii',10)
hold off
axis square
axis equal


disp([X XRec X-XRec])
