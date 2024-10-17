function plot2dplane_test
plotPoints(eye(2))

p=[-1;-2];
x0=[0.5;0.5];
q=-p'*x0;
plot2dplane(p,q);
display(q)
