%Linearly estimate homography between two sets of 2-D image points
%function H=homographyEstimateLinear(x1,x2)
%If x2 is omitted or empty, use x1(:,:,1) as points in the first view and
%x1(:,:,2:end) as points for the other view(s) (i.e., x2).
function H=homographyEstimateLinear(x1,x2)
if ~exist('x2','var') || isempty(x2)
    [x1,x2]=homographySplitData(x1);
end

NFrames=size(x2,3);
if NFrames>1
    H=zeros(3,3,NFrames);
    for iFrame=1:NFrames
        H(:,:,iFrame)=homographyEstimateLinear(x1,x2(:,:,iFrame));
    end
else
    NPoints=size(x1,2);
    x1 = homogeneous(x1,3);
    x2 = homogeneous(x2,3);

    A = zeros(2*NPoints,9);
    idxPoints=reshape(1:(2*NPoints),2,NPoints);
    for iPoint = 1:NPoints
        A(idxPoints(:,iPoint),:)=kron(x1(:,iPoint), orthComplement(x2(:,iPoint)))';
    end

    [~,~,V]=svd(A);

    H=reshape(V(:,9),3,3);
end

