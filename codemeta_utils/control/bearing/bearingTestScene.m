function sceneData=bearingTestScene(varargin)
NX=10;
L=10;
methodNoise='none';
methodTarget='random';
flagIntegerLocations=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'l'
            ivarargin=ivarargin+1;
            L=varargin{ivarargin};
        case 'nx'
            ivarargin=ivarargin+1;
            NX=varargin{ivarargin};
        case 'integerlocations'
            flagIntegerLocations=true;
        case 'sigmanoise'
            ivarargin=ivarargin+1;
            sigmaNoise=varargin{ivarargin};
        case 'noutliers'
            ivarargin=ivarargin+1;
            NOutliers=varargin{ivarargin};
        case 'methodnoise'
            ivarargin=ivarargin+1;
            methodNoise=varargin{ivarargin};
        case 'methodtarget'
            ivarargin=ivarargin+1;
            methodTarget=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

switch lower(methodTarget)
    case 'random'
        XGoal=2*L*(rand(2,1)-0.5);
    case 'southcenter'
        XGoal=[0;-L/4];
    case 'southeast'
        XGoal=[-L/3;-L/3];
    otherwise
        error('Invalid method name for target')
end

X=2*L*(rand(2,NX)-0.5);

if flagIntegerLocations
    X=round(X);
    XGoal=round(XGoal);
end

switch lower(methodNoise)
    case 'gaussian'
        XEval=X+sigmaNoise*randn(size(X));
    case 'outliers'
        XEval=[X(:,1:NX-NOutliers) 2*L*(rand(2,NOutliers)-0.5)];
    case 'none'
        XEval=X;
end

sceneData.L=2*L;
sceneData.XLandmarks=X;
sceneData.XGoal=XGoal;
sceneData.XLandmarksEval=XEval;

