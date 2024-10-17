%Returns the signed angle distance between vectors on the circle
%function a=sphere2_distSigned(y,yi)
function a=sphere2_distSigned(y,yi)
if size(yi,2)==1
    yi=permute(yi,[1 3 2]);
end

yOrth=sphere2_orth(y);
a=atan2(yOrth'*yi,y'*yi);

