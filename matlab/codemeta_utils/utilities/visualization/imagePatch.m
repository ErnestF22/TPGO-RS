%Get a patch of image around a point
%function imPatch=imagePatch(im,p,sz)
%Inputs
%   im  Input image
%   p   [2 x 1] vector with image coordinates
%   sz  [1 x 1] or [2 x 1] vector with size of the patch
%The image is padded for patches close to the edges
function imPatch=imagePatch(im,p,sz)
if numel(sz)==1
    sz=repmat(sz,2,1);
end
%transpose image coordinates into index coordinates
p=[p(2); p(1)];
halfsz=floor(sz/2);
imPadded=padarray(im,halfsz);
imPatch=imPadded(...
    p(1)+(0:sz(1)-1),...
    p(2)+(0:sz(2)-1),:);
