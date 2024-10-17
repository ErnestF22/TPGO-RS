%function d=rot_dist(R1,R2)
%Compute the geodesic distance between the rotations R1 and R2 in SO(3)

%%AUTORIGHTS%%

function d=rot_dist(R1,R2,varargin)
NR1=size(R1,3);
NR2=size(R2,3);

outputType='matrix';

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case {'vector','matrix'}
            outputType=lower(varargin{ivarargin});
        case 'outputtype'
            ivarargin=ivarargin+1;
            outputType=lower(varargin{ivarargin});
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

switch lower(outputType)
    case 'matrix'
        d=zeros(NR1,NR2);
        for r1=1:NR1
            for r2=1:NR2
                d(r1,r2)=rot_distSingle(R1(:,:,r1),R2(:,:,r2));
            end
        end
    case 'vector'
        if(NR1~=NR2)
            error('Length of inputs must be the same');
        end
        d=zeros(NR1,1);
        for r1=1:NR1
            d(r1)=rot_distSingle(R1(:,:,r1),R2(:,:,r1));
        end
        
    otherwise
        error('outputType must be ''matrix'' or ''vector''');
end

function d=rot_distSingle(R1,R2)
switch size(R1,1)
    case 2
        theta1=rot2ToAngle(R1);
        theta2=rot2ToAngle(R2);
        d=abs(modAngle(theta1-theta2));
    case 3
        R=R1'*R2;
        s1=R(6)-R(8);
        s2=R(7)-R(3);
        s3=R(2)-R(4);
        d=atan2(sqrt(s1*s1+s2*s2+s3*s3),R(1)+R(5)+R(9)-1);
    otherwise
        error('Not implemented yet')
end
%[d(r1,r2),w] = angleaxis(dc2quat(R1(:,:,r1)'*R2(:,:,r2)));
%d(r1,r2)=acos(max(min((trace(R1(:,:,r1)'*R2(:,:,r2))-1)/2,1),-1));
%d(r1,r2)=-trace(hat(logrot(R1(:,:,r1)'*R2(:,:,r2)))^2)/2;

function theta=rot2ToAngle(R)
theta=atan2(R(2,1)-R(1,2),trace(R));
