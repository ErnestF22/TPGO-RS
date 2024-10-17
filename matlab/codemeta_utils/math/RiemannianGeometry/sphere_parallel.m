%parallel transport d, a vector in the tangent space at y, along the
%geodesic given by h
function [d1]=sphere_parallel(y,h,d)
theta=norm(h);
if(theta<1e-12)
    d1=d;
else
    h=h/theta;

    dh=d'*h;
    dperp=d-h*dh;

    d1=dperp+(-y*sin(theta)+h*cos(theta))*dh;
end
