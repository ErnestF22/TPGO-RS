%Plot lines given start and end points
%function plotLines(xStart,xEnd,varargin)
%This is essentially a smart version of plot which does not require
%separating and assembling the various coordinates of the points
function plotLines(xStart,xEnd,varargin)

if isempty(varargin) || ~isStyleString(varargin{1})
    %a style has not been provided, inject default marker and size
    varargin=[{'-'} varargin{:}];
elseif styleContainsLine(varargin{1})
    %a style has been provided, but it does not contain a line style
    %so add the default one
    varargin=[{['-' varargin{1}]} varargin{2:end}];
end    

%number of input start/end points
NXStart=size(xStart,2);
NXEnd=size(xEnd,2);

%return if either is empty
if NXStart==0 || NXEnd==0
    return
end

%repeat points if start or end is a singleton
if NXStart==1 && NXEnd>1
    xStart=repmat(xStart,1,NXEnd);
    NXStart=NXEnd;
end
if NXEnd==1 && NXStart>1
    xEnd=repmat(xEnd,1,NXStart);
    NXEnd=NXStart;
end

if NXEnd~=NXStart
    error('The number of start/end points must match, or one of them must be one.')
end

sz=size(xStart);
d=sz(1);
xData=zeros(2,prod(sz(2:end)),d);
for id=1:d
    xData(:,:,id)=[squeeze(xStart(id,:,:));squeeze(xEnd(id,:,:))];
end
switch d
    case 2
        plot(xData(:,:,1),xData(:,:,2),varargin{:});
    case 3
        plot3(xData(:,:,1),xData(:,:,2),xData(:,:,3),varargin{:});
    otherwise
        error('First dimension of the data must be two or three')
end
