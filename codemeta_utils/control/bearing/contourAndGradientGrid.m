function contourAndGradientGrid(f,df,xx,yy,varargin)
flagContour=true;
flagAddSurface=false;
alphaContourAndGrad=0.3;
subFactor=10;
if ~exist('yy','var') || isempty(yy)
    yy=xx;
end
optsQuiver={};
optsContour={};
NContourLines=50;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'alpha'
            ivarargin=ivarargin+1;
            alphaContourAndGrad=varargin{ivarargin};
        case 'flagcontour'
            ivarargin=ivarargin+1;
            flagContour=varargin{ivarargin};
        case 'flagcontour'
            ivarargin=ivarargin+1;
            flagContour=varargin{ivarargin};
        case 'flagcontour'
            ivarargin=ivarargin+1;
            flagContour=varargin{ivarargin};
        case 'subfactor'
            ivarargin=ivarargin+1;
            subfactor=varargin{ivarargin};
        case 'optsquiver'
            ivarargin=ivarargin+1;
            optsQuiver=[optsQuiver varargin{ivarargin}];
        case 'optscontour'
            ivarargin=ivarargin+1;
            optsContour=[optsContour varargin{ivarargin}];
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

xxSub=xx(1:subFactor:end);
yySub=yy(1:subFactor:end);

[gridX,gridY]=meshgrid(xx,yy);
[gridXSub,gridYSub]=meshgrid(xxSub,yySub);

gridSize=size(gridX);
gridSizeSub=size(gridXSub);

C=zeros(gridSize);
D=zeros([gridSizeSub 2]);

if flagContour
    for iX=1:gridSize(1)
        for iY=1:gridSize(2)
            XEval=[gridX(iX,iY);gridY(iX,iY)];
            C(iX,iY)=f(XEval);
        end
    end
end
for iX=1:gridSizeSub(1)
    for iY=1:gridSizeSub(2)
        XEval=[gridXSub(iX,iY);gridYSub(iX,iY)];
        D(iX,iY,:)=df(XEval);
    end
end

flagHold=ishold();
if flagContour
    contour(gridX,gridY,C,NContourLines,optsContour{:})
    cm=jet(NContourLines);
    cm=1-alphaContourAndGrad+alphaContourAndGrad*cm;
    colormap(cm)
    hold on
    if flagAddSurface
        surf(gridX,gridY,C)
    end
end
quiver(gridXSub,gridYSub,-D(:,:,1),-D(:,:,2),...
    'Color',(1-alphaContourAndGrad)*[1,1,1],optsQuiver{:})
if ~flagHold
    hold off
end

axis equal
