function ny=bearingNetworkComputeRanges(x,E)
ny=squeeze(sqrt(sum((x(:,E(:,2),:)-x(:,E(:,1),:)).^2)));
