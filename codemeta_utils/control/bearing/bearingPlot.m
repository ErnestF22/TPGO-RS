%Plot bearing vectors using quiver
%function bearingPlot(X,Y,varargin)
%Inputs
%   X        base point
%   Y        bearing vectors
%   varargin all optional arguments are passed directly to quiver
function bearingPlot(X,Y,varargin)
NY=size(Y,2);
X=X*ones(1,NY);
switch size(X,1);
    case 2
        quiver(X(1,:),X(2,:),Y(1,:),Y(2,:),varargin{:})
    case 3
        quiver3(X(1,:),X(2,:),X(3,:),Y(1,:),Y(2,:),Y(3,:),varargin{:})
end        
