function [xTransformed,JRT,HRT]=projectFromG(G,X,varargin)
flagComputeJacobian=false;
flagComputeSecondJacobian=false;

if nargout>1
    flagComputeJacobian=true;
end

if nargout>2
    flagComputeSecondJacobian=true;
end

R=G2R(G);
T=G2T(G);
if ~flagComputeJacobian
    xTransformed=projectFromRT(R,T,X,varargin{:});
else
    if ~flagComputeSecondJacobian
        [xTransformed,JRT]=projectFromRT(R,T,X,varargin{:});
    else
        [xTransformed,JRT,HRT]=projectFromRT(R,T,X,varargin{:});
    end    
end
