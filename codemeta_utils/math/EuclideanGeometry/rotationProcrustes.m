%Find rotation that best aligns two sets of rotations using Frobenious norm
%function [R,R1Transformed]=rotationProcrusted(R1,R2,varargin)
%Optional arguments
%   'left'      R1=R*R2
%   'right'     R1=R2*R
function [R,R2Transformed]=rotationProcrustes(R1,R2,varargin)
NRotations=size(R1,3);
direction='left';
flagComputeTransformed=false;

if nargout>1
    flagComputeTransformed=true;
end

d=size(R1,1);

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case {'left','right'}
            direction=lower(varargin{ivarargin});
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

R12=zeros(d,d,NRotations);


switch direction
    case 'left'
        for iRotation=1:NRotations
            R12(:,:,iRotation)=R1(:,:,iRotation)*R2(:,:,iRotation)';
        end
    case 'right'
        for iRotation=1:NRotations
            R12(:,:,iRotation)=R2(:,:,iRotation)'*R1(:,:,iRotation);
        end
end

R=projectR(sum(R12,3));

if flagComputeTransformed
    R2Transformed=zeros(size(R1));
    switch direction
        case 'left'
            for iRotation=1:NRotations
                R2Transformed(:,:,iRotation)=R*R2(:,:,iRotation);
            end
        case 'right'
            for iRotation=1:NRotations
                R2Transformed(:,:,iRotation)=R2(:,:,iRotation)*R;
            end
    end
end    