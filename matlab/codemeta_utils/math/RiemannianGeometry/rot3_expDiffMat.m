%function [D,R2]=rot3_expDiffMat(R,A,varargin)
%Matrix representation of the differential of the exponential mapping as a
%Input arguments
%   R   base rotation (is a fixed parameter of the mapping exp_R)
%   A   target rotation R2 (the result of exp_R). If the 'tangent' option is
%       passed, this is a tangent vector log_R(R2)
%Optional arguments
%   'method',m
%       'closedForm'    Use method based on closed form with cotangents
%       'svd'           Use method based on SVD and closed form solution of
%                       the series defining the exponential
%Output arguments
%   D   Matrix representation of the differential of the map R2=exp_R()
%   R2  If requested, the target rotation R2
%Note: the derivation
function [D,R2]=rot3_expDiffMat(R,A,varargin)
flagAisRot=true;   %option to say that A, in fact, is Exp_R(A)
flagMethod='closedForm';    %method to use: can be 'closedForm' or 'series'

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'rot'
            flagAisRot=true;
        case 'tangent'
            flagAisRot=false;
        case 'method'
            ivarargin=ivarargin+1;
            flagMethod=varargin{ivarargin};
    end
    ivarargin=ivarargin+1;
end

if flagAisRot
    R2=A;
    A=rot_log(R,R2);
end

%pull back to the identity
Ae=R'*A;

switch lower(flagMethod)
    case 'svd'
        D=rot3_expDiffSeries(Ae);
    case 'closedform'
        [u,theta]=cnormalize(vee3(Ae));
        sincN=@(t) sinc(t/pi);
        
        D=eye(3)-theta/2*sincN(theta/2)^2*hat3(u)+(1-sincN(theta))*hat3(u)^2;
end

%compute R2 in output if requested and not already passed
if ~flagAisRot && nargout>1
    R2=rot_exp(R,A);
end
