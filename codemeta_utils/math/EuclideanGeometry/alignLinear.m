%Find the best linear transformation between two sets of vectors
%function [K,STransformed,cost]=alignLinear(S,STarget,varargin)
%Compute the matrix K such that K*S is closest to STarget
%Outputs
%   K               The transformation matrix for the alignment
%   STranformed     The transformed vectors, i.e., K*S
%   cost            The minimum of the cost, i.e., the value of ||K*S-STarget||_F^2.
%Optional arguments
%   'right'         Find K such that S*K is closest to STarget, instead of K*S
function [K,STransformed,cost]=alignLinear(S,STarget,varargin)
flagRight=false;
flagComputeSTransformed=nargout>1;
flagComputeCost=nargout>2;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'left'
            flagRight=false;
        case 'right'
            flagRight=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if flagRight
    S=S';
    STarget=STarget';
end

K=(S'\STarget')';
if flagComputeSTransformed
    STransformed=K*S;
end

if flagComputeCost
    cost=norm(STransformed-STarget,'fro')^2;
end

if flagRight
    K=K';
    if flagComputeSTransformed
        STransformed=STransformed';
    end
end
