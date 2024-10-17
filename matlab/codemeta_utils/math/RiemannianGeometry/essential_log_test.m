function essential_log_test
NPoints=100;
NGeodesics=1;
tMax=2*pi;

%resetRands()
t=linspace(0,tMax,NPoints);
err=zeros(NPoints,NGeodesics);

for iGeodesic=1:NGeodesics
    Q0=essential_randn();

    v0=[rot_randTangentNormVector(Q0(1:3,:)); rand*rot_randTangentNormVector(Q0(4:6,:))];
    v0=essential_tangentProj(Q0,v0);

    for it=1:NPoints
        Q=essential_exp(Q0,t(it)*v0);
        err(it,iGeodesic)=essential_dist(Q,essential_exp(Q0,essential_log(Q0,essential_randomVerticalMotion(Q))));
    end
end
plot(t/pi,err)
