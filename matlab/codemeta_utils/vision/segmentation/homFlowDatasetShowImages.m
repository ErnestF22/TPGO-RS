function homFlowDatasetShowImages(x,idxX)
NFrames=size(x,3);
for iFrame=1:NFrames
    subplot(1,NFrames,iFrame)
    plotGroups(x(:,:,iFrame),idxX)
    axis equal
    axis([-1 1 -1 1])
end

