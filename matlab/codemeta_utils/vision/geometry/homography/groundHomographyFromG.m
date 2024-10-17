%Compute the homography from a camera to the XY world's plane
function H=groundHomographyFromG(G1,varargin)
NPoses=size(G1,3);
flagInvertG=false;

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
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

switch methodAbsolutePoses
    case 'reference'
        flagInvertG=~flagInvertG;
    case 'pose'
        %leave flagInvertG the same
end

if flagInvertG
    G1=invg(G1);
end

%pose of a camera where the image plane corresponds to the XY plane
%in world coordinates
G0=RT2G(eye(3),[0;0;1]);

%normal of XY plane
n=[0;0;1];

%distance from plane to camera 0
d=1;

%relative rotations and translations from camera 0 to camera 1
G10=zeros(size(G1));
for iPose=1:NPoses
    G10(:,:,iPose)=G1(:,:,iPose)*invg(G0);
end

%homography
T10=G2T(G10);
n=n/d;
H=G2R(G10)+repmat(permute(T10,[1 3 2]),[1 3 1]).*repmat(n',[3 1 NPoses]);
