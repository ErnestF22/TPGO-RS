function plotPoints(x,varargin)

if isempty(varargin) || ~isStyleString(varargin{1})
    %a style has not been provided, inject default marker and size
    varargin=[{'.','MarkerSize',7} varargin{:}];
elseif styleContainsLine(varargin{1})
    %a style has been provided, but it does not contain a line style
    %so add the default one
    varargin=[['.' varargin{1}] {'MarkerSize',7} varargin{2:end}];
end    
    

d=size(x,1);
switch d
    case 2
        plot(squeeze(x(1,:,:)),squeeze(x(2,:,:)),varargin{:});
    case 3
        plot3(squeeze(x(1,:,:)),squeeze(x(2,:,:)),squeeze(x(3,:,:)),varargin{:});
    otherwise
        error('First dimension of the data must be two or three')
end




