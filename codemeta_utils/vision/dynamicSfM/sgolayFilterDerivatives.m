%Filter samples using S-G filter and compute their temporal derivatives
%function [yFilter,dyFilter,ddyFilter,g]=sgolayFilterDerivatives(dx,y,N,F)
function [yFilter,dyFilter,ddyFilter,g]=sgolayFilterDerivatives(dx,y,N,F)
% Calculate S-G coefficients
[~,g] = sgolay(N,F);   

% Resize to filter across the last dimension
sz=size(y);
y=flatten(y);

% Apply filters
yFilter=zeros(size(y));
dyFilter=zeros(size(y));
ddyFilter=zeros(size(y));
for iy=1:size(y,2)
    yFilter(:,iy)=conv(y(:,iy),g(:,1),'same');
    dyFilter(:,iy)=-conv(y(:,iy),g(:,2),'same')/dx;
    ddyFilter(:,iy)=2*conv(y(:,iy),g(:,3),'same')/(dx^2);
end

% Resize to original dimensions
yFilter=unflatten(yFilter,sz);
dyFilter=unflatten(dyFilter,sz);
ddyFilter=unflatten(ddyFilter,sz);


