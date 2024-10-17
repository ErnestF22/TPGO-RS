function [c,gradc,Hc]=essentialCost(Ri,Ti,Rj,Tj,Qij)
flagComputeGradient=false;
flagApproximatedHessian=false;
if nargout>1
    flagComputeGradient=true;
    if nargout>2
        flagApproximatedHessian=true;
    end
end

Q=essential_fromRT(Ri,Ti,Rj,Tj);
c=0.5*essential_dist(Q,Qij,'signed')^2;

if flagComputeGradient
    DQ=essential_fromRT_Diff(Ri,Ti,Rj,Tj,'Q',Q);
    gradc=-DQ'*essential_vecLog(Q,Qij,'signed');
    if flagApproximatedHessian
        DLog=essential_logDiff(Q,Qij,'components','1');
        Hc=-DQ'*DLog*DQ;
    end
end
