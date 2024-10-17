%Plot a 2-D vector field defined by a function handle
%function plotfield(f,x,y,opts)
%Plot the vector given by f([x(ix);y(iy)]) at the point [x(ix);y(iy)],
%where ix and iy go over all the elements in x and y, respectively. If y is
%omitted, y is taken to be the same as x. Options opts are directly passed
%to quiver.
function plotfield(f,x,varargin)

if ~exist('x','var')
    x=linspace(-1,1,10);
end

if ~isempty(varargin) && isnumeric(varargin{1})
    y=varargin{1};
    varargin(1)=[];
end

if ~exist('y','var')
    y=x;
end

x=shiftdim(x);
y=shiftdim(y);

if size(x,2)~=size(y,2)
    error('Dimensions of x and y should be compatible')
end

if size(x,2)==1 && size(y,2)==1
    [x,y]=meshgrid(x,y);
end

xGrid=[x(:)';y(:)'];
fxy=funEvalVec(f,xGrid);

quiver(xGrid(1,:),xGrid(2,:),fxy(1,:),fxy(2,:),varargin{:})
