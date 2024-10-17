%Draw a square patch of a 3D plane from normal
%function hOut=draw3dPlane(b,varargin)
%Input
%   b   The normal of the plane. If b is a [3 x 1] vector, the plane patch passes
%       through the origin. If b is a [4 x 1] vector, the center is
%       computed as the projection of the origin on the plane specified by
%       b in homogeneous coordinates
%Optional inputs
%   'style'             style for 'patch'
%   'side'              length of the side of the patch
%   'patchOptions',opts options to pass to patch
%   'center',p          center patch on the point p instead of the default
%
function hOut=draw3dPlane(bOrth,varargin)
flagAxesLabel=true;

NPlanes=size(bOrth,2);
if NPlanes>1
    flagHold=ishold();
    for iPlane=1:NPlanes
        draw3dPlane(bOrth(:,iPlane),varargin{:});
        hold on
    end
    if ~flagHold
        hold off
    end
else

    style='b';
    side=1;
    patchOpts={'FaceAlpha',0.5};
    o=[0;0;0];

    if size(bOrth,1)==4
        bOrth=planeNormalize(bOrth);
    end


    %optional parameters
    ivarargin=1;
    while(ivarargin<=length(varargin))
        switch(lower(varargin{ivarargin}))
            case 'style'
                ivarargin=ivarargin+1;
                style=varargin{ivarargin};
            case 'side'
                ivarargin=ivarargin+1;
                side=varargin{ivarargin};
            case 'patchopts'
                ivarargin=ivarargin+1;
                patchOpts=varargin{ivarargin};
            case 'flagaxeslabel'
                ivarargin=ivarargin+1;
                flagAxesLabel=varargin{ivarargin};
            case 'center'
                ivarargin=ivarargin+1;
                o=varargin{ivarargin};
            otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
        end
        ivarargin=ivarargin+1;
    end

    if size(bOrth,1)==4
        o=homogeneousProjectPoints(o,bOrth);
        bOrth=bOrth(1:3,:);
    end

    b1=orth([bOrth eye(3)]);
    b=b1(:,2:3);

    c1=o+side*(+b(:,1)-b(:,2));
    c2=o+side*(+b(:,1)+b(:,2));
    c3=o+side*(-b(:,1)+b(:,2));
    c4=o+side*(-b(:,1)-b(:,2));
    c=[c1 c2 c3 c4];
    hold on
    h=patch(c(1,:),c(2,:),c(3,:), style, patchOpts{:});
    if flagAxesLabel
        plot3dvect(o,bOrth,'n','b')
    else
        plot3dvect(o,bOrth,'','b')
    end
    if nargout>0
        hOut=h;
    end
end