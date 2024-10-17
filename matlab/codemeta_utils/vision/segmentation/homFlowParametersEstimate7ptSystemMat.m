%Compute the system matrix A for the LS homography flow parameter estimation
%function A=homFlowParametersEstimate7ptSystemMat(x)
function A=homFlowParametersEstimate7ptSystemMat(x)
NX=size(x,2);
A=zeros(2*NX,8);
idxA=reshape(1:2*NX,2,NX);
for iX=1:NX
    xi=x(1,iX);
    yi=x(2,iX);
    %    1  2  3  4  5  6     7     8
    Ai=[ 1 xi yi  0  0  0 xi*yi  xi^2;
         0  0  0  1 xi yi  yi^2 xi*yi;
         ];
    A(idxA(:,iX),:)=Ai;
end
