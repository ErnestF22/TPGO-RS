R0=eye(3);
R1=rot_exp(R0,rot_tangentProj(R0,0.1*randn(3)));
f1=@(x) 0.5*rot_dist(R1,x)^2;

[rot_vee(R0,-rot_log(R0,R1)) rot_vee(R0,approx_RiemGrad(rot_funs(),f1,R0))]

A=randn(3);
f2=@(x) trace(A'*x);

[rot_vee(R0,2*rot_tangentProj(R0,A)) rot_vee(R0,approx_RiemGrad(rot_funs(),f2,R0))]

A=randn(3);
f3=@(x) trace(A'*(x(:,:,1)-x(:,:,2)));
f3bb=@(R2) -trace(A'*R2);
f3b=@(R2) f3(cat(3,R1,R2));
R2=rot_exp(R0,rot_tangentProj(R0,0.1*randn(3)));
R=cat(3,R1,R2);
[rot_vee(R,2*cat(3,rot_tangentProj(R(:,:,1),A),-rot_tangentProj(R(:,:,2),A))); rot_vee(R,approx_RiemGrad(rot_funs(),f3,R))]