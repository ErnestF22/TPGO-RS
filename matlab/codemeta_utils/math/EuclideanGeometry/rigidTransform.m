%Transforms points or lines according to given rotation and translation
%function [XTransformed,JXTransformed,HXTransformed]=rigidTransform(R,T,X,varargin)
%Optional parameters
%   'references'    R,T are in the "reference" interpretation
%   'poses'         R,T are in the "pose" interpretation
%   'points'        X represents 3-D points
%   'planes'        X represents 3-D planes (must be 4-D vectors)
%Outputs
%   XTransformed    Coordinates X transformed according to R and T
%   JXTransformed   [3x6xN] matrix with Jacobians of X with respect to R and T
function [XTransformed,JXTransformed,HXTransformed]=rigidTransform(R,T,X,varargin)
flagComputeJacobian=false;
flagComputeSecondJacobian=false;
flagInvertG=false;

methodAbsolutePoses='poses';
dataType='points';
transformationDirection='wc';

if nargout>1
    flagComputeJacobian=true;
end
if nargout>2
    flagComputeSecondJacobian=true;
end

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case {'points','planes','velocities','vectors'}
            dataType=lower(varargin{ivarargin});
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=lower(varargin{ivarargin});
        case 'cw'
            transformationDirection='cw';
        case 'wc'
            transformationDirection='wc';
        case 'transformationdirection'
            ivarargin=ivarargin+1;
            transformationDirection=lower(varargin{ivarargin});
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

dimPoints=size(X,1);
NPoses=size(R,3);

if strcmp(dataType,'lines') && dimPoints~=4
    error('Lines must be represented by four-dimensional vectors')
end

switch transformationDirection
    case 'cw'
        flagInvertG=~flagInvertG;
    case 'wc'
        %leave flagInvertG the same
end

switch methodAbsolutePoses
    case 'reference'
        flagInvertG=~flagInvertG;
    case 'pose'
        %leave flagInvertG the same
end

XFlat=reshape(X,dimPoints,[]);
NPoints=size(XFlat,2);
XTransformed=zeros([size(XFlat),NPoses]);
for iPose=1:NPoses
    switch dataType
        case 'points'
            if ~flagInvertG
                XTransformed(1:3,:,iPose)=R(:,:,iPose)*XFlat(1:3,:)+T(:,iPose)*ones(1,NPoints);
            else
                XTransformed(1:3,:,iPose)=R(:,:,iPose)'*(XFlat(1:3,:)-T(:,iPose)*ones(1,NPoints));
            end
        case 'planes'
            if ~flagInvertG
                XTransformed(1:3,:,iPose)=R(:,:,iPose)*XFlat(1:3,:);
                XTransformed(4,:,iPose)=-T(:,iPose)'*XTransformed(1:3,:,iPose)+XFlat(4,:);
            else
                XTransformed(1:3,:,iPose)=R(:,:,iPose)'*XFlat(1:3,:);
                XTransformed(4,:,iPose)=T(:,iPose)'*XFlat(1:3,:)+XFlat(4,:);
            end
        case 'velocities'
            R1=R(:,:,iPose);
            T1=T(:,iPose);
            if flagInvertG
                [R1,T1]=invRT(R1,T1);
            end
            XTransformed(1:3,:,iPose)=-R1*XFlat(1:3,:);
            XTransformed(4:6,:,iPose)=hat3(XFlat(1:3,:))*R1'*T1-R1'*XFlat(4:6,:);
        case 'vectors'
            if ~flagInvertG
                XTransformed(1:3,:,iPose)=R(:,:,iPose)*XFlat(1:3,:);
            else
                XTransformed(1:3,:,iPose)=R(:,:,iPose)'*XFlat(1:3,:);
            end
        otherwise
            error('dataType is not regognized')
    end            
end

%TODO: make Jacobians for other datatype

if flagComputeJacobian
    JXTransformed=zeros(3,6,NPoints,NPoses);
    for iPose=1:NPoses
        for iX=1:NPoints
            switch methodAbsolutePoses
                case 'reference'
                    JXTransformed(:,:,iX,iPose)=[-hat3(XTransformed(:,iX,iPose)) -R(:,:,iPose)'];
                case 'pose'
                    JXTransformed(:,:,iX,iPose)=[-R(:,:,iPose)*hat(XFlat(:,iX)) eye(3)];
            end
        end
    end
end
if flagComputeSecondJacobian
    HXTransformed=zeros(6,6,3,NPoints,NPoses);
    E=eye(3);
    for iPose=1:NPoses
        for iX=1:NPoints
            XTransformedi=XTransformed(:,iX,iPose);
            Xi=XFlat(:,iX);
            for k=1:3
                ek=E(:,k);
                switch methodAbsolutePoses
                    case 'reference'
                        H1=hat(XTransformedi)*hat(ek);
                        H2=-hat(ek)*R';
                        HXTransformed(:,:,k,iX,iPose)=[H1 H2; H2' zeros(3)];
                    case 'pose'
                        H1=hat(R'*ek)*hat(Xi);
                        HXTransfxormed(:,:,k,iX,iPose)=[H1 zeros(3); zeros(3,6)];
                end
            end
        end
    end
end

XTransformed=reshape(XTransformed,[size(X) NPoses]);
