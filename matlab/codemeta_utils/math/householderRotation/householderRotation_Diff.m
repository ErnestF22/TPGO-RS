function dHVec=householderRotation_Diff(x1,x2,dx1,dx2)
if ~exist('dx2','var')
    dx2=[];
end
dx=[dx1;dx2];
DH=householderRotation_DiffMat(x1,x2);
dHVec=DH*dx;
