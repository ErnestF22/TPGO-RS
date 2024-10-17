%function R=rot_fast(vnorm,nv)
%Arguments
%   vnorm   normalized rotation axes
%   nv      rotation angles
function R=rot_fast(v,nv)
s=sin(nv);
c=cos(nv);

v1=v(1:3:end);
v2=v(2:3:end);
v3=v(3:3:end);
uc=1-c;
ucv1=uc.*v1;
ucv2=uc.*v2;
ucv1v2=ucv1.*v2;
ucv2v3=ucv2.*v3;
ucv1v3=uc.*v1.*v3;

R=[...
    c+ucv1.*v1;
    s.*v3+ucv1v2;
    -s.*v2+ucv1v3;
    -s.*v3+ucv1v2;
    c+ucv2.*v2;
    s.*v1+ucv2v3;
    s.*v2+ucv1v3;
    -s.*v1+ucv2v3;
    c+uc.*v3.^2;
    ];

