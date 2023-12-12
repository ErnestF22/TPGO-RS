%Transform R,T such that R(:,:,iNodeFix)=eye(d) and T(:,iNodeFix)=zeros(d,1)
%function [R,T]=RTFix(R,T,iNodeFix)
function [R,T]=RTFix(R,T,iNodeFix,varargin)
methodAbsolutePoses='poses';

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
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

RFix=R(:,:,iNodeFix);
TFix=T(:,iNodeFix);

switch methodAbsolutePoses
    case 'reference'
        R=multiprod(RFix',R);
        T=RFix'*(T-TFix);
    case 'pose'
        error('Not implemented')
    otherwise
        error('methodAbsolutePoses not recognized')
end