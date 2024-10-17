%function essential_kernelQuickShiftTreeSpy(t)
%Visualizes the adjacency matrix of a tree produced by essential_kernelQuickShiftTree
function essential_kernelQuickShiftTreeSpy(t)
N=length(t);
tVis=zeros(N);
for iN=1:N
    if t(iN)~=0
        tVis(iN,t(iN))=1;
    end
end
imagesc(tVis)