function [Ti,A,b]=sfm_rawAverageTranslationsDirect(Ri,Tij,E,varargin)
NEdges=size(E,1);
flagCenter=true;
flagRescale=false;

methodAbsolutePoses='reference';

%Get the number of rotations from maximum number in edges
%Can be changed by passing corresponding argument
NTranslations=max(E(:));

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'ntranslations'
            ivarargin=ivarargin+1;
            NTranslations=varargin{ivarargin};
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=lower(varargin{ivarargin});
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end


idxI=zeros(18,NEdges);
idxJ=zeros(18,NEdges);
idxS=zeros(18,NEdges);
b=zeros(3,NEdges);
idxTranslation=reshape(1:3*NTranslations,[3 NTranslations]);
I3=eye(3);

for iEdge=1:NEdges
    iTranslation=E(iEdge,1);
    jTranslation=E(iEdge,2);
    
    idxIdx=(3*iEdge-2:3*iEdge)';
    
    iIdxTranslation=idxTranslation(:,iTranslation)';
    jIdxTranslation=idxTranslation(:,jTranslation)';
    
    idxI(:,iEdge)=repmat(idxIdx,6,1);
    idxJ(:,iEdge)=reshape(repmat([iIdxTranslation jIdxTranslation],3,1),18,1);

    switch methodAbsolutePoses
        case 'reference'
            RiTranslation=Ri(:,:,iTranslation);
            idxS(:,iEdge)=reshape([-I3 I3],18,1);
            b(:,iEdge)=RiTranslation*Tij(:,iEdge);
        case 'pose'
            RiTranslation=Ri(:,:,iTranslation);
            RjTranslation=Ri(:,:,jTranslation);
            idxS(:,iEdge)=reshape([I3 RiTranslation'*RjTranslation],18,1);
            b(:,iEdge)=Tij(:,iEdge);
        otherwise
            error('Invalid method for absolute poses')
    end
end

A=sparse(idxI,idxJ,idxS);
b=b(:);

%A is rank deficient by construction
%We need to use the SVD
[U,S]=svd(full(A'*A));

TiVec=U(:,1:end-3)...
    *(diag(diag(S(1:end-3,1:end-3)).^-1) ...
    *(U(:,1:end-3)'*(A'*b)));
Ti=reshape(TiVec,3,[]);

if flagCenter
    Ti=Ti-mean(Ti,2)*ones(1,NTranslations);
end
if flagRescale
    Ti=Ti/norm(Ti(:));
end

