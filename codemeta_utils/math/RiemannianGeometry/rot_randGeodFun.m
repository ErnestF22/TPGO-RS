%function [Rt,vt,R0,v0,vVec]=rot_randGeodFun(A,varargin)
%Same as rot_geodFun, but the rotations and tangent vector for the geodesic
%are generated at random. 
%If A is omitted, generates a geodesic in SO(3)
%If A is a scalar, generates a geodesic in SO(n), n=A
%If A is a 2-D matrix, generate a geodesic starting from R=A
%If A is a 3-D matrix, generate multiple geodesics, each  starting from
%   R(:,:,iN), iN=1,...,N
%Warning: The function does not check if A is valid rotation
%
%Optional key/value arguments
%   'perp', w   where w is a tangent vector at the origin of the
%               geodesic. In this case the generated vector will be
%               orthogonal to w.
%   'speed',s   generate a geodesic with speed s (by default, s=1)
%   'randSpeed' random geodesic speed in [0,1], same as 'speed',rand
%
%See also rot_geodFun
function [Rt,dRt,R0,dR0,vVec,ddRt,dvVec]=rot_randGeodFun(A,varargin)
flagAWasNotRotation=false;  %remember if the passed argument A was a matrix
s=1;                        %speed of the geodesic
if ~exist('A','var') || isempty(A)
    A=3;
end
if isscalar(A)
    R=orth(randn(A));
    R=R/det(R);
    flagAWasNotRotation=true;
else
    R=A;
end

perpTangent=[];

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'perp'
            if flagAWasNotRotation
                error(['You must specify the origin of the geodesic if you want to also specify an orthognal direction'])
            end
            ivarargin=ivarargin+1;
            perpTangent=varargin{ivarargin};
        case 'speed'
            ivarargin=ivarargin+1;
            s=varargin{ivarargin};
        case 'randspeed'
            s=rand;
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if isempty(perpTangent)
    v=rot_randTangentNormVector(R);
else
    v=rot_randTangentPerpNormVector(R,perpTangent);
end
NR=size(A,3);
if length(s)==1 && NR>1
    s=s*ones(1,NR);
end
s=permute(s,[1 3 2]);
[Rt,dRt,R0,dR0,vVec,ddRt,dvVec]=rot_geodFun(R,v,'speed',s);
