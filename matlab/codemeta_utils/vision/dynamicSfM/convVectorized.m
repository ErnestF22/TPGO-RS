function [yg,idxValid]=convVectorized(y,g,varargin)
flagComputeIdx=nargout>1;

sz=size(y);
y=flatten(y);
[ly,Ny]=size(y);
ygTest=conv(y(:,1),g,varargin{:});
lyg=length(ygTest);
yg=zeros(lyg,Ny);
for iy=1:Ny
    yg(:,iy)=conv(y(:,iy),g,varargin{:});
end
sz(end)=lyg;
yg=unflatten(yg,sz);
if flagComputeIdx
    idxValid=convValidIdx(ly,length(g),varargin{:});
end

