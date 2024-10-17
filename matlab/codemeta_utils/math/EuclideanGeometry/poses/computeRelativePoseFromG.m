%Compute relative pose from two absolute poses
%function G12=computeRelativePoseFromG(G1,G2,X2,varargin)
%The transformation G12 transforms the data in the reference frame of G2 to
%the reference frame of G1
%Optional inputs
%   '21'    give the transformation from the frame of G2 to the one of G1
%           (default)
%   '12'    give the transformation from the frame of G1 to the one of G2
function [G12]=computeRelativePoseFromG(G1,G2,varargin)
methodAbsolutePoses='reference';
flagInvertG=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case '12'
            flagInvertG=true;
        case '21'
            flagInvertG=false;
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

N1=size(G1,3);
N2=size(G2,3);
G12=zeros(4,4,N1,N2);
for iN1=1:N1
    for iN2=1:N2
        switch methodAbsolutePoses
            case 'reference'
                G12(:,:,iN1,iN2)=invg(G1(:,:,iN1))*G2(:,:,iN2);
            case 'pose'
                G12(:,:,iN1,iN2)=G1(:,:,iN1)*invg(G2(:,:,iN2));
            otherwise
                error('Absolute pose method specification not valid')
        end
    end
end
G12=squeeze(G12);
if flagInvertG
    G12=invg(G12);
end
