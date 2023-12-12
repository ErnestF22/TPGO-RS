function POCPreRigidRotationContractionStabilityConditionQuadratic
syms m11 m12 m22 lg k c b
ai=-2*k*m12+2*m11*b;
ci=-k*m22+m11-c*m12+2*m12*b;
di=m12+(b-c)*m22;

q=ai*di-ci^2;
[c,v]=coeffs(q,{m11,m12});

Aq=[c(1) c(2)/2; c(2)/2 c(4)];
bq=[c(3); c(5)];

m=[m11;m12];
disp('The following is zero if Aq and bq are correct')
disp(simplify(m.'*Aq*m+bq.'*m+c(6)-q))
m0=-inv(Aq)*bq;
disp(m0)
disp(c(6))
