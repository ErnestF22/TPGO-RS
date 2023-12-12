function POCDerivativeGeodesicQuadraticSpeed
resetRands();
s=@(t) t^2/2;
ds=@(t) t;
dds=@(t) 1;
R0=rot_randn();
wHat=rot_hat(R0,randn(3,1));
R=@(t) rot_exp(R0,s(t)*wHat);
dR=@(t) ds(t)*R(t)*R0'*wHat;
ddR=@(t) dds(t)*R(t)*R0'*wHat+ds(t)^2*R(t)*R0'*wHat^2;

funCheckDer(R,dR,'angle')
