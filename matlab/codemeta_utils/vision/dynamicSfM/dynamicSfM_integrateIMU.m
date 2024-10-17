
function [R,T]=dynamicSfM_integrateIMU(t,wIMU,alphaIMU,R0,T0,v0)
g=gravityVector();

R=rotDyn_integrateVelocity(t,wIMU,R0,'left');
a=multiprodMatVec(invR(R),alphaIMU)-repmat(g,1,size(alphaIMU,2));
T=realDyn_integrateAcceleration(t,a,T0,v0);
