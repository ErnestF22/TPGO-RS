function l2d=project3DLineTo2DLine(R,T,l3d)
%Note: l3d must is assumed to be a [4 x 2] matrix
%pose interpretation
NLines=size(l3d,3);
l2d=zeros(3,NLines);
for iLine=1:NLines
    A=([[R'; T'] l3d(:,:,iLine)]);
    c=null(A);
    l2d(:,iLine)=c(1:3);
end
