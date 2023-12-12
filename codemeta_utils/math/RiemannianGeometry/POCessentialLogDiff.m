function POCessentialLogDiff

Q1=essential_eye();
v1=essential_randTangentNormVector(Q1);
Q2=essential_exp(Q1,0.1*v1);

[Q2t,vt,Q20,v0,vVec]=essential_randGeodFun(Q2);

u=@(t) essential_vee(Q1,essential_log(Q1,Q2t(t)));

u1=@(t) [eye(3) zeros(3)]*u(t);
u2=@(t) [zeros(3) eye(3)]*u(t);

Q11=essential_getR1(Q1);
Q12=essential_getR2(Q1);
Q21=@(t) essential_getR1(Q2t(t));
Q22=@(t) essential_getR2(Q2t(t));

Q2tClosest=@(t) essential_closestRepresentative(Q1,Q2t(t));
Q21Closest=@(t) essential_getR1(Q2tClosest(t));
Q22Closest=@(t) essential_getR2(Q2tClosest(t));

t=rand; disp([Q21(t) Q11*rot(u1(t))])

Du=@(t) blkdiag(...
    rot3_logDiff(Q11,Q21Closest(t)),...
    rot3_logDiff(Q12,Q22Closest(t)));

du=@(t) Du(t)*vVec;

t=linspace(-pi/8,pi/8,30);
%plotfun(@(t) essential_distMinAngle(Q2t(0),Q2t(t)),t)
check_der(u,du,t)
