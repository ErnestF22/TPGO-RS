%Compute relative pose from absolute poses
%function [R12,T12]=computeRelativePoseFromRT(R1,T1,R2,T2,varargin)
%Returns the rigid change of reference (R12,T12) from camera 2 to camera 1.
function [R12,T12]=computeRelativePoseFromRT(R1,T1,R2,T2,varargin)
methodAbsolutePoses='reference';
flagInvertDirection=false;

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
        case 'invertdirection'
            flagInvertDirection=true;
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if flagInvertDirection
    tempR=R1;
    tempT=T1;
    R1=R2;
    T1=T2;
    R2=tempR;
    T2=tempT;
end

switch methodAbsolutePoses
    case 'reference'
        R12=R1'*R2;
        T12=R1'*(T2-T1);
    case 'pose'
        R12=R1*R2';
        T12=T1-R12*T2;
    otherwise
        error('Absolute pose method specification not valid')
end
