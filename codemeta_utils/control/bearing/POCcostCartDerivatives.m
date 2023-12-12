function POCcostCartDerivatives
Xg=randn(3,1);
Xe0=randn(3,1);
dXe=cnormalize(randn(3,1));
Xi=randn(3,1);
Xe=@(t) Xe0+t*dXe;

Yei=@(t) cnormalize(Xi-Xe(t));
nYei=@(t) norm(Xi-Xe(t));
Ygi=cnormalize(Xi-Xg);
ci=@(t) min(1,max(-1,Ygi'*Yei(t)));
ai=@(t) acos(ci(t));
Jei=@(t) (eye(3)-Yei(t)*Yei(t)')/nYei(t);
Hi=@(t) -(Yei(t)*Ygi'*Jei(t)+Jei(t)*Ygi*Yei(t)'+ci(t)*Jei(t))/nYei(t);
dci=@(t) -Ygi'*Jei(t)*dXe;
ddci=@(t) dXe'*Hi(t)*dXe;
%check_der(ci,dci)
%check_der(dci,ddci)

dai=@(t) -dci(t)/sqrt(1-ci(t)^2);
ddai=@(t) -1/sqrt(1-ci(t)^2)^3*ci(t)*dci(t)^2-1/sqrt(1-ci(t)^2)*ddci(t);
%check_der(ai,dai)
%check_der(dai,ddai)

gradPhi=@(t) ai(t)/sqrt(1-ci(t)^2)*Jei(t)*Ygi;
HessPhi=@(t) ((1/(1-ci(t)^2)-ai(t)*ci(t)/sqrt(1-ci(t)^2)^3)*Jei(t)*Ygi*Ygi'*Jei(t)...
    -ai(t)/sqrt(1-ci(t)^2)*Hi(t));
phi=@(t) 0.5*ai(t)^2;
dphi=@(t) dXe'*gradPhi(t);
%ddphi=@(t) dXe'*Jei(t)*Ygi*Ygi'*Jei(t)*dXe/(1-ci(t)^2)...
%    -ai(t)*(1/sqrt(1-ci(t)^2)^3*ci(t)*dXe'*Jei(t)*Ygi*Ygi'*Jei(t)*dXe...
%    +1/sqrt(1-ci(t)^2)*dXe'*Hi(t)*dXe);
% ddphi=@(t) dXe'*(Jei(t)*Ygi*Ygi'*Jei(t)/(1-ci(t)^2)...
%     -ai(t)*ci(t)/sqrt(1-ci(t)^2)^3*Jei(t)*Ygi*Ygi'*Jei(t)...
%     -ai(t)/sqrt(1-ci(t)^2)*Hi(t))*dXe;
ddphi=@(t) dXe'*HessPhi(t)*dXe;

%check_der(phi,dphi)
check_der(dphi,ddphi)
