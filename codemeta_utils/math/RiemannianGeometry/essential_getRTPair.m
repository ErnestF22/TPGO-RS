%Get pair of cameras corresponding to a point in the QREM
%function [R1,T1,R2,T2]=essential_getRTPair(Q,varargin)
%Inputs
%   Q     a point in the QREM
%Optional inputs
%   'references'    (R,T)'s are in the "reference" interpretation
%   'poses'         (R,T)'s are in the "pose" interpretation
function [R1,T1,R2,T2]=essential_getRTPair(Q,varargin)
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
R1=essential_getR1(Q);
T1=zeros(3,1);
R2=essential_getR2(Q);
T2=[0;0;1];

if strcmpi(methodAbsolutePoses,'pose')
    [R1,T1]=invRT(R1,T1);
    [R2,T2]=invRT(R2,T2);
end
