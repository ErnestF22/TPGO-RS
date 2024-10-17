%Apply a 2D affine transformation to feature
%function x=sfm_rawTransformFeature(K,x,varargin)
function x=sfm_rawTransformFeature(K,x,varargin)
flagInvert=false;       %invert the transformation

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'invert'
            flagInvert=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

switch numel(K)
    case 1
        %K is a scalar
        if ~flagInvert
            x=x*K;
        else
            x=x/K;
        end
    case 4
        %K is an homogeneous transformation
        if ~flagInvert
            x=K*x;
        else
            x=K\x;
        end
    case 6
        %K is an affine transformation
        Nx=size(x,2);
        if ~flagInvert
            x=K(1:2,1:2)*x+K(1:2,3)*ones(1,Nx);
        else
            x=K(1:2,1:2)\(x-K(1:2,3)*ones(1,Nx));
        end
    case 9
        %K is a projective transformation
        Nx=size(x,2);
        if ~flagInvert
            x1=K*[x(1:2,:); ones(1,Nx)];
        else
            x1=K\[x(1:2,:); ones(1,Nx)];
        end
        x=x1(1:2,:)./([1;1]*x1(3,:));
end
