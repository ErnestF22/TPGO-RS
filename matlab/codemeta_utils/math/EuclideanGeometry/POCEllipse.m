function POCEllipse
U=rot_randn(eye(2));
S=diag([1 2]);
V=rot_randn(eye(2))';
M=U*S*V';
a=U*S;
t=linspace(0,2*pi,100);
x=[cos(t); sin(t)];
Mx=M*x;
plot(x(1,:),x(2,:))
hold on
plot(Mx(1,:),Mx(2,:),'r')
plot(a(1,:),a(2,:),'go')
hold off
