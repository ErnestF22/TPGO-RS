function POCWLSasGaussNewton
global xt;
global Asqrt;
global dx;

preA=randn(3);
A=preA'*preA;
Asqrt=sqrtPSDMatrix(A);
[xt,~,~,dx]=real_randGeodFun(randn(3,1));

f1=@(t) 0.5*xt(t)'*A*xt(t);
df1=@(t) xt(t)'*A*dx;
ddf1=@(t) dx'*A*dx;

%check_der(df1,ddf1)

%funCompare(f1,@f2)
funCompare(df1,@df2)
%funCompare(ddf1,@ddf2)



function f=f2(t)
global xt;
global Asqrt;
f=0;
E=eye(3);
for k=1:3
    f=f+0.5*(E(:,k)'*Asqrt*xt(t))^2;
end

function f=df2(t)
global xt
global Asqrt
global dx
f=0;
E=eye(3);
for k=1:3
    %f=f+(E(:,k)'*Asqrt*xt(t))*E(:,k)'*Asqrt*dx;
    f=f+xt(t)'*Asqrt*E(:,k)*E(:,k)'*Asqrt*dx;
end

function f=ddf2(t)
global xt
global Asqrt
global dx
f=0;
E=eye(3);
for k=1:3
    f=f+(E(:,k)'*Asqrt*dx)*E(:,k)'*Asqrt*dx;
end


function funCompare(f1,f2)
t=linspace(0,1,100);
plotfun(f1,t)
hold on
plotfun(f2,t,'rx')
hold off
