function bearingCostDisplay(sceneData,funs,varargin)
NGrid=91;
optsX={'MarkerSize',7};
optsGoal={'MarkerFaceColor','g'};
optsQuiver={};
optsContour={};

flagLandmarksEval=false;
flagContourAndGrad=true;
flagContour=true;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'ngrid'
            ivarargin=ivarargin+1;
            NGrid=varargin{ivarargin};
        case 'flagcontourandgrad'
            ivarargin=ivarargin+1;
            flagContourAndGrad=varargin{ivarargin};
        case 'flagcontour'
            ivarargin=ivarargin+1;
            flagContour=varargin{ivarargin};
        case 'optsx'
            ivarargin=ivarargin+1;
            optsX=[optsX varargin{ivarargin}];
        case 'optsgoal'
            ivarargin=ivarargin+1;
            optsGoal=[optsGoal varargin{ivarargin}];
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

L=sceneData.L;
XGoal=sceneData.XGoal;
if isfield(sceneData,'XLandmarksEval')
    flagLandmarksEval=true;
    XEval=sceneData.XLandmarksEval;
    XOriginal=sceneData.XLandmarks;
else
    XOriginal=sceneData.XLandmarks;
end
    

YGoal=bearingCompute(XGoal,XOriginal);

xGrid=linspace(-L/2,L/2,NGrid);

if flagContourAndGrad
    if ~flagLandmarksEval
        contourAndGradientGrid(@(x) cost(x,XOriginal,YGoal,funs), @(x) gradCost(x,XOriginal,YGoal,funs),xGrid,[],...
            'flagContour',flagContour,'optsContour',optsContour,'optsQuiver',optsQuiver);
    else
        contourAndGradientGrid(@(x) cost(x,XEval,YGoal,funs), @(x) gradCost(x,XEval,YGoal,funs),xGrid,[],...
            'flagContour',flagContour,'optsContour',optsContour,'optsQuiver',optsQuiver);
    end
end

hold on
plot(XOriginal(1,:),XOriginal(2,:),'b*',optsX{:})
plot(XGoal(1),XGoal(2),'go',optsX{:},optsGoal{:})
if flagLandmarksEval
    plot(XEval(1,:),XEval(2,:),'bo',optsX{:})
    plot([XEval(1,:); XOriginal(1,:)],[XEval(2,:); XOriginal(2,:)],'b:')
end
hold off
axis([-L/2 L/2 -L/2 L/2])
axis square
axis equal

function c=cost(XEval,X,YGoal,funs)
[YEval,nYEval]=bearingCompute(XEval,X);
c=bearingCostGeneral(YEval,YGoal,nYEval,funs);

function gradc=gradCost(XEval,X,YGoal,funs)
YEval=bearingCompute(XEval,X);
gradc=bearingCostGeneral_gradient(YEval,YGoal,funs);
