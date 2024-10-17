%Project 3-D points and Jacobians of R, T
%[xTransformed,JRT,HRT]=projectFromRT(R,T,X,varargin)
%If NX=size(X,2) and NP=size(R,3)=size(T,2), then
%   xTransformed    [2 x NX x NP] matrix with the projected image points
%   JRT             [2 x 6 x NX x NP] matrix with the Jacobian of the
%                   projections w.r.t. the rotations and translations
%   HRT             [6 x 6 x 2 x NX x NP] matrix with the second Jacobian
%                   of the projections w.r.t. the rotations and translations
%If NX==1, the the matrices are squeezed (e.g., xTransformed becomes [2 x NP]
%
%TODO: fix second Jacobian
function [xTransformed,JRT,HRT]=projectFromRT(R,T,X,varargin)
flagComputeJacobian=false;
flagComputeSecondJacobian=false;
flagReshapedX=false;

if nargout>1
    flagComputeJacobian=true;
end

if nargout>2
    flagComputeSecondJacobian=true;
end

NX3=size(X,3);
if NX3>1
    NX2=size(X,2);
    X=X(:,:);           %flatten X, we will reshape the result at the end
    flagReshapedX=true;
end

NPoses=size(R,3);
NX=size(X,2);
if NPoses>1
    %do recursive call for multiple poses
    xTransformed=zeros(2,NX,NPoses);
    if flagComputeJacobian
        JRT=zeros(2,6,NX,NPoses);
        if flagComputeSecondJacobian
            HRT=zeros(6,6,2,NX,NPoses);
        end
    end
    for iPose=1:NPoses
        if ~flagComputeJacobian
            xTransformed(:,:,iPose)=projectFromRT(R(:,:,iPose),T(:,iPose),X,varargin{:});
        else
            if ~flagComputeSecondJacobian
                [xTransformed(:,:,iPose),JRT(:,:,:,iPose)]=projectFromRT(R(:,:,iPose),T(:,iPose),X,varargin{:});
            else
                [xTransformed(:,:,iPose),JRT(:,:,:,iPose),HRT(:,:,:,:,iPose)]=projectFromRT(R(:,:,iPose),T(:,iPose),X,varargin{:});
            end    
        end
    end
    xTransformed=squeeze(xTransformed);
    if flagComputeJacobian
        JRT=squeeze(JRT);
        if flagComputeSecondJacobian
            HRT=squeeze(HRT);
        end
    end
else
    methodAbsolutePoses='pose';
    
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

    if ~flagComputeJacobian
        XTransformed=rigidTransform(R,T,X,'methodAbsolutePoses',methodAbsolutePoses);
    else
        if ~flagComputeSecondJacobian
            [XTransformed,JXTransformed]=rigidTransform(R,T,X,'methodAbsolutePoses',methodAbsolutePoses);
        else
            [XTransformed,JXTransformed,HXTransformed]=rigidTransform(R,T,X,'methodAbsolutePoses',methodAbsolutePoses);
        end
    end

    xTransformed=XTransformed(1:2,:)./repmat(XTransformed(3,:),2,1);
    if flagComputeJacobian
        JRT=zeros(2,6,NX);
        for iX=1:NX
            JRT(:,:,iX)=([eye(2) -xTransformed(:,iX)]*JXTransformed(:,:,iX))./XTransformed(3,iX);
        end
    end
    if flagComputeSecondJacobian
        HRT=zeros(6,6,2,NX);
        E=eye(3);
        e3=E(:,3);
        for iX=1:NX
            Xi=XTransformed(:,iX);
            Ji=JXTransformed(:,:,iX);
            Hi=HXTransformed(:,:,:,iX);
            e3Xi=XTransformed(3,iX);

            f1=1/e3Xi^2;
            M=(e3'*Xi*eye(3)-Xi*e3');
            for k=1:2
                ek=E(:,k);
                H1=((-2/e3Xi)*Ji'*e3*ek'*M*Ji);
                H2=Ji'*(e3*ek'-ek*e3')*Ji;
                H3=zeros(6,6);
                for l=1:3
                    el=E(:,l);
                    H3=H3+ek'*M*el*Hi(:,:,l);
                end
                HRT(:,:,k,iX)=f1*(H1+H2+H3);
            end
    %         for k=1:2
    %             ek=E(:,k);
    %             eSym=(e3'*Xi)*ek-(ek'*Xi)*e3;
    %             H1=zeros(6,6);
    %             for k2=1:3
    %                 H1=H1+eSym(k2)*HXTransformed(:,:,k2,iX);
    %             end
    %             HRT(:,:,k,iX)=(H1*e3Xi-2*Ji'*eSym*e3'*Ji)/e3Xi^3;
    %         end
        end
    end
end

%   xTransformed    [2 x NX x NP] matrix with the projected image points
%   JRT             [2 x 6 x NX x NP] matrix with the Jacobian of the
%                   projections w.r.t. the rotations and translations
%   HRT             [6 x 6 x 2 x NX x NP] matrix with the second Jacobian
%                   of the projections w.r.t. the rotations and translations
if flagReshapedX
    xTransformed=reshape(xTransformed,[2 NX2 NX3 NPoses]);
    if flagComputeJacobian
        JRT=reshape(JRT,[2 6 NX2 NX3 NPoses]);
        if flagComputeSecondJacobian
            HRT=reshape(JRT,[6 6 2 NX NPoses]);
        end
    end
end

        
end