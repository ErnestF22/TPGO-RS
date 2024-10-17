function [l,dl,Hl]=logLikelihood(Ri,Rj,Ti,Tj,Rij,Tij,Gamma,EType)
flagComputeGradient=false;
flagApproximatedHessian=false;
if nargout>1
    flagComputeGradient=true;
    if nargout>2
        flagApproximatedHessian=true;
    end
end
if ~exist('EType','var')
    EType=1;
end
switch EType
    case 1
        vRij=logrot(Rij'*Ri'*Rj);
        vTij=Tij-Ri'*(Tj-Ti);
    case 2
        vRij=logrot(Rij'*Ri*Rj');
        vTij=Tij-(Ti-Ri*Rj'*Tj);
end
v=[vRij;vTij];
l=(v'*Gamma*v)/2;

if flagComputeGradient
    Dij=rot3_logDiff(Rij,vRij,'tangentVec');
    GammaRR=Gamma(1:3,1:3);
    GammaRT=Gamma(1:3,4:6);
    GammaTT=Gamma(4:6,4:6);
    
    switch EType
        case 1
            hij=hat(Ri'*(Tj-Ti));

            gR1=GammaRR*vRij;
            gR2=GammaRT*vTij;
            gR=Dij'*(gR1+gR2);

            gT1=GammaTT'*vTij;
            gT2=GammaRT'*vRij;
            gT=gT1+gT2;
            gTb=Ri*gT;

            dl=[
                -Ri'*Rj*gR+hij*gT gR; ...
                gTb -gTb
            ];
        case 2
            gR_R=Rj'*Dij'*(GammaRR*vRij+GammaRT*vTij);
            gTCommon=GammaTT*vTij+GammaRT'*vRij;
            gTCommon2=Ri'*gTCommon;
            gT_Ti=-gTCommon;
            gT_Tj=Rj*gTCommon2;
            gT_R=hat(Rj'*Tj)*gTCommon2;
            gR=gR_R+gT_R;
            dl=[gR -gR; gT_Ti gT_Tj];
    end
    if flagApproximatedHessian
        if l==0
            Hl=zeros(12);
        else
            Hl=(dl(:)*dl(:)')/((2*l)^2);
        end
    end
end