function POCHessianDistanceRot3

[R1,dR1,~,~,v1]=rot_randGeodFun(eye(3));
[R2,dR2,~,~,v2]=rot_randGeodFun(eye(3),'speed',rand);

%funCheckDer(@(t) costAndDer(R1(t),R2(t),v1,v2),'function',linspace(0,10))
%funCheckDer(@(t) derAndDder(R1(t),R2(t),v1,v2),'function',linspace(0,10))

funPlot(@(t) eig(hess(R1(t),R2(t))))

function H=hess(R1,R2)
D=rot3_logDiff(R1,R2);
symD=multisym(D);
H=[symD -D; -D' symD];

function [d,dd]=costAndDer(R1,R2,v1,v2)
d=0.5*rot_dist(R1,R2)^2;
dd=-logrot(R1'*R2)'*v1-logrot(R2'*R1)'*v2;

function [dd,ddd]=derAndDder(R1,R2,v1,v2)
dd=-logrot(R1'*R2)'*v1-logrot(R2'*R1)'*v2;
ddd=[v1;v2]'*hess(R1,R2)*[v1;v2];



