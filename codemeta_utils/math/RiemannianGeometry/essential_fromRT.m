%Compute a point in the QREM corresponding to a pair of cameras
%function Q=essential_fromRT(R1,T1,R2,T2,varargin)
%Inputs
%   R1,T1,R2,T2     camera poses
%Optional inputs
%   'references'    (R,T)'s are in the "reference" interpretation
%   'poses'         (R,T)'s are in the "pose" interpretation

function Q=essential_fromRT(R1,T1,R2,T2,varargin)
methodAbsolutePoses='reference';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=lower(varargin{ivarargin});
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

N=size(R1,3);
if ~exist('R2','var') || isempty(R2)
    R2=R1;
    T2=T1;
    R1=repmat(eye(3,3),[1 1 N]);
    T1=zeros(3,N);
end

if strcmpi(methodAbsolutePoses,'pose')
    [R1,T1]=invRT(R1,T1);
    [R2,T2]=invRT(R2,T2);
end


v=T2-T1;
R0=householderRotation(v,3);

Q=zeros(6,3,N);
for iN=1:N
    Q(1:3,1:3,iN)=R0(:,:,iN)*R1(:,:,iN);
    Q(4:6,1:3,iN)=R0(:,:,iN)*R2(:,:,iN);
end
