%function x2=homographyMap(H,x1)
%Maps points in view one to points in view two as x2~H*x1
%Optional inputs
%   'inverse'   use inv(H) instead of H
%   'line'      use the homography for lines instead of coordinates (i.e.,
%               H' instead of H
function x2=homographyMap(H,x1,varargin)

N=size(H,3);
if N>1
    x2=zeros([size(x1) N]);
    for iN=1:N
        x2(:,:,iN)=homographyMap(H(:,:,iN),x1,varargin{:});
    end
else
    flagInverse=false;
    flagLine=false;

    %optional parameters
    ivarargin=1;
    while(ivarargin<=length(varargin))
        switch(lower(varargin{ivarargin}))
            case 'inverse'
                flagInverse=true;
            case 'line'
                flagLine=true;
            otherwise
                    error(['Argument ' varargin{ivarargin} ' not valid!'])
        end
        ivarargin=ivarargin+1;
    end
    x1Hom=homogeneous(x1,3);

    if ~flagLine
        if ~flagInverse
            x2Hom=H*x1Hom;
        else
            x2Hom=H\x1Hom;
        end
    else
        if ~flagInverse
            x2Hom=H'\x1Hom;
        else
            x2Hom=H'*x1Hom;
        end
    end

    x2=homogeneous(x2Hom,2);
end
