%function plotSphereMap(f,varargin)
%Plots a spherical coodinate grid as deformed by a map f
%Inputs
%   f   a function from [3 x 1] vectors to [3 x 1] vectors. If omitted or empty,
%       the identity map
%Optional arguments
%   'NRadii',NRadii             number of steps on the radial coordinate
%   'NLatitutes',NLatitudes     number of steps on the latitude coordinate
%   'NLongitudes',NLongitudes   number of steps on the longitudes
%                               coordinate
%   'stepRadii',stepRadii       step size for the radial coordinate
function plotSphereMap(f,varargin)

if ~exist('f','var') || isempty(f)
    f=@(x) x;
end

[p,NRadii,NLatitudes,NLongitudes]=sphereGrid(varargin{:});
p=permute(p,[2 3 4 1]);
p=evalfunVec(f,p);

cmap=winter(NRadii);
for iRadii=1:NRadii
    plotArgs={'b-','color',cmap(iRadii,:)};
    if iRadii~=NRadii
        pPlot=reshape(p(iRadii:iRadii+1,:,:,:),2,[],3);
        plot3(pPlot(:,:,1),pPlot(:,:,2),pPlot(:,:,3),plotArgs{:})
    end
    hold on
    for iLatitudes=1:NLatitudes-1
        pPlot=reshape(permute(p(iRadii,iLatitudes:iLatitudes+1,:,:),[2 1 3 4]),2,[],3);
        plot3(pPlot(:,:,1),pPlot(:,:,2),pPlot(:,:,3),plotArgs{:})
    end
    for iLongitudes=1:NLongitudes
        idxLongitudes=mod((iLongitudes:iLongitudes+1)-1,NLongitudes)+1;
        pPlot=reshape(permute(p(iRadii,:,idxLongitudes,:),[3 1 2 4]),2,[],3);
        plot3(pPlot(:,:,1),pPlot(:,:,2),pPlot(:,:,3),plotArgs{:})
    end    
end
hold off
axis equal
