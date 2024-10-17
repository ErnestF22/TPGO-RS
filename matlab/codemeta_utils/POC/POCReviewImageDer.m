function POCReviewImageDer
x=sym('x',[3 1]);
vr=sym('vr',[3 1]);
w=sym('w',[3 1]);

dx=vr+cross(w,x);

y=[x(1)/x(3);x(2)/x(3);1/x(3)];

dydx=[diff(y,x(1)) diff(y,x(2)) diff(y,x(3))];

dy=dydx*dx;

yp=sym('yp',[3 1]);
xp=[yp(1)/yp(3); yp(2)/yp(3); 1/yp(3)];

dyp=simplify(subs(dy,x,xp));

collect(expand(dyp),[w;vr])

%simplify(dyp(1)-(-yp(1)*yp(2)*w(1)+(1+yp(1)^2)*w(2))+yp(2)*w(3))

%simplify(dyp(2)-(+yp(1)*yp(2)*w(2)-(1+yp(2)^2)*w(1))-yp(1)*w(3))









