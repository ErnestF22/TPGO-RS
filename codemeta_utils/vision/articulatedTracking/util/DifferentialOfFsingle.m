function DF = DifferentialOfFsingle(x,v)

tol = 10^(-5);
if abs(norm(x)-1)>tol
    error('x should be a unit vector...')
end



nv = max(norm(v),eps);
u = v/nv;


[F1,F2] = F(x,v);

I = eye(3);
F11 = cos(nv)*I;
F21 = -nv*sin(nv)*I;
F12 = -sin(nv)*x*u' + cos(nv)*(u*u') + (sin(nv)/nv)*(I - (u*u'));
F22 = -( (sin(nv)/nv) + cos(nv))*x*v' + cos(nv)*I - nv*sin(nv)*(u*u');

DF =[F11 F12;F21 F22];

DF = [eye(3) zeros(3);F1*F2' eye(3)]*DF*[eye(3) zeros(3); -x*v' eye(3)];



end

