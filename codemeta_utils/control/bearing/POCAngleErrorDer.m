function POCAngleErrorDer
resetRands();
theta10=modAngle(randn);
dtheta1=rand;
theta20=modAngle(randn);
dtheta2=rand;

theta1=@(t) theta10+t*dtheta1;
theta2=@(t) theta20+t*dtheta2;

f=@(t) modAngle(theta1(t)-theta2(t));
df=@(t) dtheta1-dtheta2;

check_der(f,df,linspace(-20,20))

