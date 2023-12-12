%function [R,T]=planeToRT(NVec,varargin)
%Return a reference frame with the origin in the plane and the z axis
%aligned to the normal of the plane.
function [RPlane,TPlane]=planeToRT(NVec,varargin)
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
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NPlanes=size(NVec,2);
if NPlanes>1
    RPlane=zeros(3,3,NPlanes);
    TPlane=zeros(3,NPlanes);
    for iPlane=1:NPlanes
        [RPlane(:,:,iPlane),TPlane(:,iPlane)]=planeToRT(NVec(:,iPlane),varargin{:});
    end
else
    [N,d]=planeNVecToNd(NVec);
    RPlane=fliplr(orthCompleteBasis(N));
    %flip first column to maintain positive determinant
    RPlane(:,1)=-RPlane(:,1);
    TPlane=d*N;

    switch methodAbsolutePoses
        case 'reference'
            %nothing to do
        case 'pose'
            [RPlane,TPlane]=invRT(RPlane,TPlane);
        otherwise
            error('methodAbsolutePoses invalid')
    end
end
