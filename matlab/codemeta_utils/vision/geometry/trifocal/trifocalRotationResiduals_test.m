function trifocalRotationResiduals_test

Nn=2;
n0=squeeze(sphere_randn([],[],Nn));
n1=squeeze(sphere_randn([],[],Nn));
n2=squeeze(sphere_randn([],[],Nn));

[R1,~,~,~,v1Vec]=rot_randGeodFun([],'speed',1);
[R2,~,~,~,v2Vec]=rot_randGeodFun([],'speed',rand);

f=@(t) funDer(R1(t),R2(t),n0,n1,n2,v1Vec,v2Vec);
g=@(t) gradDerGrad(R1(t),R2(t),n0,n1,n2,v1Vec,v2Vec);

%check_der(f)
check_der(g)

function [e,de]=funDer(R1,R2,n0,n1,n2,v1Vec,v2Vec)
v12Vec=[v1Vec;v2Vec];
[e,ge]=trifocalRotationResiduals(R1,R2,n0,n1,n2);
de=v12Vec'*ge;

function [ge,dge]=gradDerGrad(R1,R2,n0,n1,n2,v1Vec,v2Vec)
v12Vec=[v1Vec;v2Vec];
[~,ge,Dge]=trifocalRotationResiduals(R1,R2,n0,n1,n2);
Nn=size(ge,2);
dge=zeros(size(ge));
for in=1:Nn
    dge(:,in)=Dge(:,:,in)*v12Vec;
end
