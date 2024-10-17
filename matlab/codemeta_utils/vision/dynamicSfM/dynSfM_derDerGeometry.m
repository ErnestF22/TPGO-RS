%Second derivative of geometry (3-D points) expressed in the coordinates of a moving camera
function ddXb=dynSfM_derDerGeometry(Xb,wb,dwb,nub,alphab)
NX=size(Xb,2);
hwb=hat3(wb);
hwbSq=multiprod(hwb,hwb);
hdwb=hat3(dwb);
ddXb=multiprod(hwbSq-hdwb,Xb)+repmat(permute(2*multiprodMatVec(hwb,nub)-alphab,[1 3 2]),1,NX,1);

