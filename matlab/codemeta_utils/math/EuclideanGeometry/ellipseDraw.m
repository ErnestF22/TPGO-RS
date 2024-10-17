%Draw a ellipse given focus and major radius a
function ellipseDraw(focus,a)
l = diff(focus,[],2);
theta = atan(l(2)/l(1));
center = sum(focus,2)/2;
c = norm(diff(focus,[],2))/2;
b = sqrt(a^2-c^2);
t = linspace(0,pi*2);
x =@(t) a*cos(t);
y =@(t) b*sin(t);
nx =@(t) x(t)*cos(theta)-y(t)*sin(theta)+center(1);
ny =@(t) x(t)*sin(theta)+y(t)*cos(theta)+center(2);
plot(nx(t),ny(t),'-.r')


