%Projects stacked approximate rotations to rotation matrices
%function Ri=sfm_rawRotationsVecProject(RiVec,varargin)
function Ri=sfm_rawRotationsVecProject(RiVec,varargin)
methodAbsolutePoses='reference';

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
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NRotations=size(RiVec,1)/3;
Ri=reshape(RiVec',[3 3 NRotations]);

signs=zeros(NRotations,1);
for iRotation=1:NRotations
    signs(iRotation)=sign(det(Ri(:,:,iRotation)));
end
Ri=Ri*median(signs);

Ri=projectR(Ri);

switch methodAbsolutePoses
    case 'reference'
        %nothing to do
    case 'pose'
        Ri=invR(Ri);
    otherwise
        error('Invalid method for absolute poses');
end
