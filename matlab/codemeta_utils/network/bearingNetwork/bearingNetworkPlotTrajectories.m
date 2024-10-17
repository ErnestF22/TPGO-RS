function bearingNetworkPlotTrajectories(x,style,varargin)
if ~exist('style','var')
    style='-';
end
flagPlotInitialPositions=true;
color=0.75*ones(1,3);

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'color'
            ivarargin=ivarargin+1;
            color=varargin{ivarargin};
        case 'flagplotinitialpositions'
            ivarargin=ivarargin+1;
            flagPlotInitialPositions=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end


switch size(x,1)
    case {2,4}
        d=2;
    case {3,6}
        d=3;
    otherwise
        error('Dimension of the space not recognized')
end

if flagPlotInitialPositions
    xInit=x(1:d,:,1);
    plotPoints(xInit,'rx')
end
flagHold=ishold();
hold on
N=size(x,2);
if size(color,1)==1
    color=ones(N,1)*color;
end

switch d
    case 2
        for iN=1:N
            plot(squeeze(x(1,iN,:))',squeeze(x(2,iN,:))',style,'color',color(iN,:))
        end
    case 3
        for iN=1:N
            plot3(squeeze(x(1,iN,:))',squeeze(x(2,iN,:))',squeeze(x(3,iN,:))',...
                style,'color',color(iN,:))
        end
end

if ~flagHold
    hold off
end
