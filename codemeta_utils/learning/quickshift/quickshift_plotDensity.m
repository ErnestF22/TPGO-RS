%Plot the density used by QuickShift
%function quickshift_plotDensity(X,phi,varargin)
%Input parameters
%   X       [2 x NPoints] matrix with the datapoints
%   phi     Function handle to the density kernel
%Optional parameters
%   'optsDensity',opts  opts is a cell array of options to pass to
%       quickshift_density
%   'limits', ax    ax is a [1 x 4] array with the limits of where the
%       density should be plotted. If omitted, these are determined by the
%       values of X
function quickshift_plotDensity(X,phi,varargin)
optsDensity={};
minWidth=1;
flagAutoLimits=true;


%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'optsdensity'
            ivarargin=ivarargin+1;
            optsDensity=varargin{ivarargin};
        case 'limits'
            ivarargin=ivarargin+1;
            ax=varargin{ivarargin};
            XMin=ax([1,3]);
            XMax=ax([2,4]);
            flagAutoLimits=false;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if flagAutoLimits
    %compute grid size
    XMin=min(X,[],2);
    XMax=max(X,[],2);
    XWidth=max(XMax-XMin,minWidth);
    XMin=XMin-XWidth/4;
    XMax=XMax+XMax/4;
end

%compute grid points
NGrid=100;
[xGrid,yGrid]=meshgrid(linspace(XMin(1),XMax(1),NGrid),linspace(XMin(2),XMax(2),NGrid));
XGrid=[xGrid(:)';yGrid(:)'];
%compute distances
D=sqrt(euclideanDistMatrix(X,[X XGrid]));
P=quickshift_density(phi,D,optsDensity{:});
NPoints=size(X,2);
PGrid=P(NPoints+1:end);
P=P(1:NPoints);
plot3(X(1,:),X(2,:),P,'ro','MarkerFaceColor','r');
hold on
surf(reshape(XGrid(1,:),NGrid,NGrid),...
    reshape(XGrid(2,:),NGrid,NGrid),...
    reshape(PGrid,NGrid,NGrid),'EdgeColor','None');
hold off
%axis equal
