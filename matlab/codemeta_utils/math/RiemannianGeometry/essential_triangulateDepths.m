function l=essential_triangulateDepths(Q,x)
x=homogeneous(x,3);
R1=essential_getR1(Q);
R2=essential_getR2(Q);

Nx=size(x,2);
l=zeros(2,Nx);
for ix=1:Nx
    A=[R1*x(:,ix,1) -R2*x(:,ix,2)];
    l(:,ix)=(A'*A)\(A'*[0;0;1]);
end
