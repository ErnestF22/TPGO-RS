function POCTransformationRelativeAngle

%methodTrajectories='cross';
methodTrajectories='randlimited';
%methodTrajectories='oppositeHorizontal';
%methodAngle='bearingVerticalProjection';
methodAngle='signed';
TFinal=5;


switch lower(methodAngle)
    case 'cosine'
        fAngle=@(x1,x2) cnormalize(x1)'*cnormalize(x2);
    case 'unsigned'
        fAngle=@(x1,x2) vctAngle(x1,x2);
    case 'signed'
        fAngle=@(x1,x2) atan2(x2(2),x2(1))-atan2(x1(2),x1(1));
    case 'bearingverticalprojection'
        fAngle=@(x1,x2) orthComplement(x1)'*cnormalize(x2);
end

switch lower(methodTrajectories)
    case 'oppositevertical'
        x1=real_geodFun([1;1],[0;-1]);
        x2=real_geodFun([1;-1],[0;1]);
    case 'oppositehorizontal'
        x1=real_geodFun([1;1],[-1;0]);
        x2=real_geodFun([-0.5;-1],[1;0]);
    case 'cross'
        x1=real_geodFun([0.3;0.2],cnormalize([1;0.3]));
        x2=real_geodFun([0;1],cnormalize([1;-0.3]));
    case 'rand'
        %resetRands(2)
        x1=real_randGeodFun(rand(2,1));
        x2=real_randGeodFun(rand(2,1));
    case 'randlimited'
        %resetRands(6)
        flagGenerate=true;
        while flagGenerate
            x1=real_randGeodFun(rand(2,1));
            x2=real_randGeodFun(rand(2,1));
            a0=fAngle(x1(0),x2(0));
            if abs(a0)>deg2rad(45) && abs(a0)<deg2rad(90)
                flagGenerate=false;
            end
        end
        display(rad2deg(a0))
        
end



c=@(t) fAngle(x1(t),x2(t));

A=@(x1,x2) householderRotation(x1,1);
At=@(t) A(x1(t),x2(t));

x1Transformed=@(t) At(t)*x1(t);
x2Transformed=@(t) At(t)*x2(t);
%cTransformed=@(t) vctAngle(x1Transformed(t),x2Transformed(t));
cTransformed=@(t) -fAngle(x1Transformed(t),x2Transformed(t));

t=linspace(0,TFinal,100);
subplot(2,2,1)
plotTrajectories(x1,x2,t)
subplot(2,2,2)
plotTrajectories(x1Transformed,x2Transformed,t)
subplot(2,2,3)
funPlot(c,t,'bx')
hold on
funPlot(cTransformed,t,'r')
hold off

c0=c(0);

e=@(t) abs(c(t)-c0);

subplot(2,2,4)
funPlot(e,t)
ylabel('Angular error [rad]')
%set(gca,'YLim',[0,2*pi])
grid on

function plotTrajectories(x1,x2,t)
funPlotCurve(x1,t,'bx')
hold on
funPlotCurve(x2,t,'r-')
plot(0,0,'ko')
hold off
axis equal

function funPlotCurve(x,t,varargin)
if ~exist('t','var') || isempty(t)
    t=linspace(0,1,100);
end
xt=funEval(x,t);
plotPoints(xt,varargin{:})

