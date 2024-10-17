%function Y=lie_randn(funs,y,v,N)
%Generate random tangent vectors at y from a isotropic zero mean Gaussian
%distribution with variance v and return the corresponding element in the
%manifold
%If N>1, generate N rotations in this way and return them in a 3-D array 
%by stacking all the matrices along the 3rd dimension
%Default values if argument are omitted or empty
%   y   funs.eye
%   v   100
%   N   1
function Y=lie_randn(funs,y,v,N)
if ~exist('y','var') || isempty(y)
    y=funs.eye();
end
Ny=size(y,3);

if ~exist('v','var') || isempty(v)
    v=100;
end

if ~exist('N','var') || isempty(N)
    N=Ny;
elseif N~=Ny && Ny~=1
    error('MATLAB:argument','Number of elements in y should be one or N');
end

y0=y(:,:,1);
if N>1
    Y=zeros([size(y,1) size(y,2) N]);
    for iN=1:N
        if Ny>1
            y0=y(:,:,iN);
        end
        Y(:,:,iN)=lie_randn(funs,y0,v,1);
    end
else
    Y=funs.exp(y,v*randn*funs.randTangentNormVector(y));
end
