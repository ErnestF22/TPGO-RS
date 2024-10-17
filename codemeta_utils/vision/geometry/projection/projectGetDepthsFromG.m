%Get the depths (lambda) of points seen by the given camera
%function lambda=projectGetDepthsFromG(G,X,varargin)
function [lambda,JRT]=projectGetDepthsFromG(G,X,varargin)
flagComputeJacobian=nargout>1;
[R,T]=G2RT(G);
if flagComputeJacobian
    [lambda,JRT]=projectGetDepthsFromRT(R,T,X,varargin{:});
else
    lambda=projectGetDepthsFromRT(R,T,X,varargin{:});
end    

