function showHingeEnsambleResponse(w,v,a,b,varargin)
L=1;
NxxPlot=100;
flagDisplayResponse=false;
flagDisplayContour=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'l'
            ivarargin=ivarargin+1;
            L=varargin{ivarargin};
        case 'n'
            ivarargin=ivarargin+1;
            NxxPlot=varargin{ivarargin};
        case 'displayresponse'
            flagDisplayResponse=true;
        case 'displaycontour'
            flagDisplayContour=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

xxPlot=linspace(-L,L,NxxPlot);
[xPlot,yPlot]=meshgrid(xxPlot,xxPlot);
XPlot=reshape(permute(cat(3,xPlot,yPlot),[3 2 1]),2,[]);

r=hingeEnsambleResponse(XPlot,w,v,a,b);
r=reshape(r,NxxPlot,NxxPlot);

imagesc(xxPlot,xxPlot,r)
colormap gray

if flagDisplayResponse
    disp(r)
end

if flagDisplayContour
    %hold on
    contour(xxPlot,xxPlot,r,[0 0])
    %hold off
end


