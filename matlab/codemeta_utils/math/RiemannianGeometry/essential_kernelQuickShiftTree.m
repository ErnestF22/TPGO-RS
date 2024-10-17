%function t=essential_kernelQuickShiftTree(d)
%Given a distance matrix d, and a kernel estimate f, create the QuickShift tree
%Inputs
%   d   [N x N] matrix distances, might contain NaNs
%   f   [N x 1] vector with values of the kernel estimate
%Outputs
%   t   [N x 1] vector containing the index of the parent. t==0 for the
%       root
function t=essential_kernelQuickShiftTree(d,f)
N=length(f);
t=zeros(N,1);
for iN=1:N
    idxValid=find(~isnan(d(:,iN)) & f>f(iN));
    if ~isempty(idxValid)
        [~,idxMin]=min(d(iN,idxValid));
        t(iN)=idxValid(idxMin);
    end
end
