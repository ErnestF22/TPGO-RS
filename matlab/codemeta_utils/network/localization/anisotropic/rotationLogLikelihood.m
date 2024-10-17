function [l,dl,Hl]=rotationLogLikelihood(Ri,Rj,Rij,Gamma)
flagComputeGradient=false;
flagApproximatedHessian=false;
if nargout>1
    flagComputeGradient=true;
    if nargout>2
        flagApproximatedHessian=true;
    end
end

v=logrot(Rij'*Ri'*Rj);
Gammav=Gamma*v;
l=(v'*Gammav)/2;

if flagComputeGradient
    D2=rot3_logDiff(Rij,v,'tangentVec')';
    D1=-Ri'*Rj*D2;
    dl=reshape([D1;D2]*Gammav,3,2);
    if flagApproximatedHessian
        if l==0
            Hl=zeros(6);
        else
            Hl=[D1;D2]*Gamma*[D1;D2]';
        end
    end
end
