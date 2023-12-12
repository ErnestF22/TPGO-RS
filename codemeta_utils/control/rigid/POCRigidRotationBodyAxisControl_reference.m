function [er,der,dder]=POCRigidRotationBodyAxisControl_reference(A,T)
if ~exist('T','var')
    T=1;
end

%generate reference normal vector (er) and its derivatives
RT=rot_randn();
a=@(t) A*sin(2*pi/T*t);
da=@(t) A*2*pi/T*cos(2*pi/T*t);
dda=@(t) -A*4*pi^2/T^2*sin(2*pi/T*t);
%check_der(a,da)
%check_der(da,dda)

xr=@(t) cos(a(t));
dxr=@(t) -sin(a(t))*da(t);
ddxr=@(t) -cos(a(t))*da(t)^2-sin(a(t))*dda(t);
yr=@(t) sin(a(t));
dyr=@(t) cos(a(t))*da(t);
ddyr=@(t) -sin(a(t))*da(t)^2+cos(a(t))*dda(t);

% check_der(xr,dxr)
% check_der(dxr,ddxr)
% check_der(yr,dyr)
% check_der(dyr,ddyr)

er=@(t) RT*[xr(t);yr(t);0];
der=@(t) RT*[dxr(t);dyr(t);0];
dder=@(t) RT*[ddxr(t);ddyr(t);0];
%plotPoints(evalfun(er))
% check_der(er,der)
% check_der(der,dder)
