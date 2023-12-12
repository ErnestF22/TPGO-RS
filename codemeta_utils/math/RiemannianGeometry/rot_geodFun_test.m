function rot_geodFun_test

[Rt,dRt,R0,dR0,vVec,ddRt,dvVec]=rot_geodFun([],[]);

figure(1)
funCheckDer(Rt,dRt,'angle')
figure(2)
funCompare(vVec, @(t) rot_vee(Rt(t),dRt(t)))
figure(3)
funCheckDer(dRt,ddRt,'angle')
figure(4)
funCheckDer(vVec,dvVec,'angle')

