function [Ti,A]=sfm_rawAverageTranslationsFromRotationsAndPoints(Ri,x1,x2,E,varargin)
NEdges=size(E,1);
flagCenter=true;
flagRescale=true;
flagSignFromDepths=true;
%flagL1=true;

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

switch methodAbsolutePoses
    case 'reference'
        %nothing to do
    case 'pose'
        Ri=invR(Ri);
    otherwise
        error('Invalid method for absolute poses')
end

NPoints=cellfun(@(x) size(x,2),x1);
NPointsIdx=[1;cumsum(NPoints)+1];
NPointsTot=sum(NPoints);

idxTranslation=reshape(1:3*NTranslations,[3 NTranslations]);

idxJ=zeros(6,NPointsTot);
idxI=zeros(6,NPointsTot);
idxS=zeros(6,NPointsTot);

for iEdge=1:NEdges
    iTranslation=E(iEdge,1);
    jTranslation=E(iEdge,2);
    
    idxIdx=NPointsIdx(iEdge):(NPointsIdx(iEdge+1)-1);
    
    RiTranslation=Ri(:,:,iTranslation);
    RjTranslation=Ri(:,:,jTranslation);
    
    coeffs=cross(...
        RiTranslation*homogeneous(x1{iEdge},3),...
        RjTranslation*homogeneous(x2{iEdge},3)...
        );
            
    idxJ(:,idxIdx)=repmat(idxIdx,6,1);
    idxI(:,idxIdx)=repmat([idxTranslation(:,iTranslation);idxTranslation(:,jTranslation)],1,NPoints(iEdge));
    idxS(:,idxIdx)=[coeffs;-coeffs];
end

A=sparse(idxJ,idxI,idxS);

% if ~flagL1
[U,S,V]=svd(full(A'*A));
Ti=reshape(V(:,end-3),[3 NTranslations]);
% else
%     cvx_begin
%         variable Ti(3,NTranslations)
%         minimize(norm(A*Ti(:),1))
%         subject to
%             sum(Ti,2)==[0;0;0]
%     cvx_end
% end
    
if flagCenter
    Ti=Ti-mean(Ti,2)*ones(1,NTranslations);
end
if flagRescale
    Ti=Ti/norm(Ti(:));
end
if flagSignFromDepths
    l=zeros(2,NPointsTot);
    for iEdge=1:NEdges
        iTranslation=E(iEdge,1);
        jTranslation=E(iEdge,2);

        RiTranslation=Ri(:,:,iTranslation);
        RjTranslation=Ri(:,:,jTranslation);
        
        TiTranslation=Ti(:,iTranslation);
        TjTranslation=Ti(:,jTranslation);

        for iPoint=1:NPoints(iEdge)
            ALambda=[ RiTranslation*homogeneous(x1{iEdge}(:,iPoint),3) ...
               -RjTranslation*homogeneous(x2{iEdge}(:,iPoint),3)];
            
            bLambda=TjTranslation-TiTranslation;
            l(:,NPointsIdx(iEdge)+iPoint-1)=ALambda\bLambda;
        end
    end
    Ti=Ti*median(sign(l(:)));
end


switch methodAbsolutePoses
    case 'reference'
        %nothing to do
    case 'pose'
        [~,Ti]=invRT(Ri,Ti);
    otherwise
        error('Invalid method for absolute poses')
end


