function POClogLikelihoodDerivativeType2

[Ri,~,~,~,vRi]=rot_randGeodFun([],'speed',rand);
[Rj,~,~,~,vRj]=rot_randGeodFun([],'speed',rand);
[Ti,~,~,vTi]=real_randGeodFun(randn(3,1),'speed',rand);
[Tj,~,~,vTj]=real_randGeodFun(randn(3,1),'speed',rand);

preGamma=randn(6);
Gamma=preGamma'*preGamma;
Rij=(Ri(0)*Rj(0)')';
Tij=Ti(0)-Ri(0)*Rj(0)'*Tj(0);

%RRes=@(t) residualRAndDer(Ri(t),Rj(t),Rij,vRi,vRj);
%check_der(RRes)
% TRes=@(t) residualTAndDer(Ri(t),Ti(t),Rj(t),Tj(t),Tij,vRi,vTi,vRj,vTj);
% check_der(TRes)
% vRes=@(t) logResidualRAndDer(Ri(t),Rj(t),Rij,vRi,vRj);
% check_der(vRes)
l=@(t) costAndDer(Ri(t),Ti(t),Rj(t),Tj(t),Rij,Tij,Gamma,vRi,vTi,vRj,vTj);
check_der(l)



function [R,dR]=residualRAndDer(Ri,Rj,Rij,vRi,vRj)
R=Rij*Ri*Rj';
dR=Rij*Ri*Rj'*hat(Rj*(vRi-vRj));

function [vR,dvR]=logResidualRAndDer(Ri,Rj,Rij,vRi,vRj)
vR=logrot(Rij'*Ri*Rj');
Dij=rot3_logDiff(Rij,vR,'tangentVec');
dvR=Dij*(Rj*(vRi-vRj));

function [vT,dvT]=residualTAndDer(Ri,Ti,Rj,Tj,Tij,vRi,vTi,vRj,vTj)
vT=Tij-(Ti-Ri*Rj'*Tj);
dvT=-vTi-Ri*hat(Rj'*Tj)*(vRi-vRj)+Ri*Rj'*vTj;

function [l,dl]=costAndDer(Ri,Ti,Rj,Tj,Rij,Tij,Gamma,vRi,vTi,vRj,vTj)
vRij=logrot(Rij'*Ri*Rj');
vTij=Tij-(Ti-Ri*Rj'*Tj);
Dij=rot3_logDiff(Rij,vRij,'tangentVec');
dvR=Dij*Rj*(vRi-vRj);
dvT=-vTi-Ri*hat(Rj'*Tj)*(vRi-vRj)+Ri*Rj'*vTj;
v=[vRij;vTij];
l=(v'*Gamma*v)/2;

GammaRR=Gamma(1:3,1:3);
GammaRT=Gamma(1:3,4:6);
GammaTT=Gamma(4:6,4:6);

% dl=(vR'*GammaRR+vT'*GammaRT')*dvR+(vT'*GammaTT+vR'*GammaRT)*dvT;
% dl=(vRij'*GammaRR+vTij'*GammaRT')*Dij*Rj*(vRi-vRj)...
%     +(vTij'*GammaTT+vRij'*GammaRT)*(-vTi-Ri*hat(Rj'*Tj)*(vRi-vRj)+Ri*Rj'*vTj);
gR_R=Rj'*Dij'*(GammaRR*vRij+GammaRT*vTij);
gTCommon=GammaTT*vTij+GammaRT'*vRij;
gTCommon2=Ri'*gTCommon;
gT_Ti=-gTCommon;
gT_Tj=Rj*gTCommon2;
gT_R=hat(Rj'*Tj)*gTCommon2;
gR=gR_R+gT_R;
% dl=gR'*(vRi-vRj)...
%     +gT_Ti'*vTi+gT_R'*(vRi-vRj)+gT_Tj'*vTj;
gradl=[gR -gR; gT_Ti gT_Tj];
dl=trace([vRi vRj; vTi vTj]'*gradl);
