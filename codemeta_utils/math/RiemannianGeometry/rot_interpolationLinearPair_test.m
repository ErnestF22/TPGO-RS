function rot_interpolationLinearPair_test
switch 3
    case 1
        R1=rot_randn();
        R2=rot_randn();
    case 2
        R1=eye(3);
        R2=rot([0;0;1]*pi/4);
    case 3
        R1=rot([0;0;1]*pi/4);
        R2=R1*rot([0;1;0]*pi/4);
end

v=rot_vee(R1,rot_log(R1,R2));
[u,d]=cnormalize(v);

t0=rand;
t1=t0+rand;

t=linspace(t0,t1,100);

RInterp=rot_interpolationLinearPair([t0 t1],cat(3,R1,R2),t);

vInterp=rot_vee(R1,rot_log(R1,RInterp));
[uInterp,dInterp]=cnormalize(vInterp);

eu=diag(euclideanDistMatrix(u*dInterp,vInterp));

figure(1)
subplot(3,1,1)
plot(eu)
subplot(3,1,2)
plot(dInterp/d)
subplot(3,1,3)
plotRotationTrajectory(RInterp);
