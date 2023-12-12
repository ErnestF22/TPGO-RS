function POCPreRigidRotationContractionStabilityConditions
%note k and lg can be absorbed together into k
%either m11 or m22 can be fixed (because the conditions are homogeneous in M
%%
syms m11 m12 m22 lg k c b
ai=-2*k*m12+2*m11*b;
ci=-k*m22+m11-c*m12+2*m12*b;
di=m12+(b-c)*m22;

constraints={...
    ai<0,...
    di<0,...
    ai*di-ci^2>0,...
    m11>0,...
    m22>0,...
    m11*m22-m12^2>0};
%%    
NGrid=10;
figure(1)
xGrid=linspace(0,3,NGrid);
yGrid=xGrid;
substitutions={m22,1,k,1,c,1.8,b,1};
Nconstraints=length(constraints);
fplot=@(cn) funImageConstraints(cn,m11,m12,substitutions,...
    'xGrid',xGrid,...
    'yGrid',yGrid);
 

q=ai*di-ci^2;
[cq,vq]=coeffs(q,{m11,m12});

Aq=[cq(1) cq(2)/2; cq(2)/2 cq(4)];
bq=[cq(3); cq(5)];
m=[m11;m12];
disp(simplify(m.'*Aq*m+bq.'*m+cq(6)-q))

Aq=double(subs(Aq,substitutions(1:2:end),substitutions(2:2:end)));
bq=double(subs(bq,substitutions(1:2:end),substitutions(2:2:end)));

m0=-Aq\bq/2;

qsub=subs(q,substitutions(1:2:end),substitutions(2:2:end));
qfun=@(v) double(subs(qsub,{'m11' 'm12'},{v(1) v(2)}));
%funImage(xGrid,yGrid,qfun);

display(m0)
display(qfun(m0))

for ic=1:Nconstraints
    subplot(1,Nconstraints+1,ic)
    fplot(constraints(ic))
    hold on
    plotPoints(m0,'x')
    hold off
    axis equal
end
subplot(1,Nconstraints+1,Nconstraints+1)
fplot(constraints)
axis equal
hold on
plotPoints(m0,'x')
hold off

fplot(constraints)
axis equal

% disp('Cond1')
% disp(simplify(ai<0))
% disp('Cond2')
% disp(di<0)
% disp('Cond3')
% disp(collect(simplify(expand(ai*di-ci>0)),[m12 m22]))


