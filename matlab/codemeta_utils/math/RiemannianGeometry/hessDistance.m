function h=hessDistance(wi,wj,Sij)

if(exist('Sij','var')==0)
    Sij=eye(3);
end
bracket=@(A,B) A*B-B*A;
phiij=norm(logrot(Sij));
if(phiij>2*eps)
    c1=(1-phiij*cos(phiij)/sin(phiij))/(2*sin(phiij)^2);
    c2=phiij/sin(phiij);
else
    c1=0;
    c2=1;
end
hwi=hat(wi);
hwj=hat(wj);

h=c1*trace((hwi-hwj)*Sij)^2+c2*trace(((hwi-hwj)*(hwi-hwj)'+bracket(hwi,hwj))*Sij);
