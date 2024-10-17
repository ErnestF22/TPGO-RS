%Show and animate graph with a linear deformation
function graphAnimate(E,x,dx,R,dR)
if exist('R','var') && ~isempty(R)
    flagDrawLocalAxes=true;
    RAnimated=R;
else
    flagDrawLocalAxes=false;
    R=[];
    RAnimated=[];
    
end
if exist('dR','var') && ~isempty(dR)
    flagAnimateLocalAxes=true;
else
    flagAnimateLocalAxes=false;
end

timesDuration=5;
timesFPS=15;
deformationAmplitude=1;
deformationPeriod=1;

timesLength=timesDuration*timesFPS;
times=linspace(0,timesDuration,timesLength);
deformationCoeffs=deformationAmplitude*sin(2*pi*times/deformationPeriod);
graphShow(E,x,R)
ax=axis;
for iTimes=1:timesLength
    xAnimated=x+deformationCoeffs(iTimes)*dx;
    if flagAnimateLocalAxes
        RAnimated=rot_exp(R,rot_hat(R,deformationCoeffs(iTimes)*dR));
    end
    graphShow(E,xAnimated,RAnimated)
    axis(ax)
    pause(1/timesFPS);
end

