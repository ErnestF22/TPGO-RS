%Plot scales as circles centered at the corresponding points
function quickshift_plotScales(X,scales,varargin)
NPoints=size(X,2);
for iPoint=1:NPoints
    plotCircle(X(1,iPoint),X(2,iPoint),scales(iPoint),varargin{:});
end
