function POCStudyControlCart
g0=randn(2,1);
dg=randn(2,1);
g=@(t) g0+t*dg;

thetad=@(t) atan2([0 1]*g(t),[1 0]*g(t));
dthetad=@(t) (g(t)'*[0 1; -1 0]*dg)/(g(t)'*g(t));

check_der(thetad, dthetad)