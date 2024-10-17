%Derivative of geometry (3-D points) expressed in the coordinates of a moving camera
function dXb=dynSfM_derGeometry(Xb,wb,nub)
NX=size(Xb,2);
hwb=hat3(wb);
dXb=-multiprod(hwb,Xb)-repmat(permute(nub,[1 3 2]),1,NX,1);
