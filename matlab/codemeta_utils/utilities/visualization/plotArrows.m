%Plot arrows given start and end points
%function plotArrows(xStart,xEnd,varargin)
%This is essentially a smart version of plot which does not require
%separating and assembling the various coordinates of the points
function plotArrows(xStart,xEnd,varargin)

if isempty(varargin) || ~isStyleString(varargin{1})
    %a style has not been provided, inject default marker and size
    varargin=[{'-'} varargin{:}];
elseif styleContainsLine(varargin{1})
    %a style has been provided, but it does not contain a line style
    %so add the default one
    varargin=[{['-' varargin{1}]} varargin{2:end}];
end    

if isempty(xStart)
    return
end

%automatically expand one of the arguments if necessary
NS=size(xStart,2);
NE=size(xEnd,2);
if NS~=NE
    if NS==1
        xStart=repmat(xStart,1,NE);
    end
    if NE==1
        xEnd=repmat(xEnd,1,NS);
    end
end
        

sz=size(xStart);
d=sz(1);
switch d
    case 2
        quiver(xStart(1,:),xStart(2,:),...
            xEnd(1,:)-xStart(1,:),xEnd(2,:)-xStart(2,:),0,varargin{:});
    case 3
        quiver3(xStart(1,:),xStart(2,:),xStart(3,:),...
            xEnd(1,:)-xStart(1,:),xEnd(2,:)-xStart(2,:),xEnd(3,:)-xStart(3,:),0,varargin{:});
    otherwise
        error('First dimension of the data must be two or three')
end
