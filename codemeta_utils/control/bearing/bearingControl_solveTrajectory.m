%Simulates the trajectory of an agent using bearing control
%function [x,t]=bearingControl_solveTrajectory(sceneData,funs,varargin)
%Inputs
%   sceneData   Struct with fields
%       XLandmarksEval  location of landmarks used for the computation of
%                       the bearings
%       XLandmarksEval  location of landmarks used for the computation of
%                       the goal bearings
%       XGoal           location of goal
%       L               length of a side of the scene
%   funs        Reshaping functions generated with bearingCostFunctions
%Optional arguments
%   'x0',x0         Initial coordinates
%   'tFinal',tF     Simulated time for which to run the simulation for
%   'model',str     Which kinematic model to use. Can be:
%       'direct'    single integrator
%       'unicycle'  a 2-D kinematic cart (unicycle model). In this case the
%           state includes a third entry with the angle
%
%Note: it would probably be better to rewrite this in Simulink       
function [x,t]=bearingControl_solveTrajectory(sceneData,funs,varargin)
%short-hand notation to scene data
XLandmarks=sceneData.XLandmarks;
XLandmarksEval=sceneData.XLandmarksEval;
XGoal=sceneData.XGoal;
if ~isa(XGoal, 'function_handle')
    XGoal=@(t) XGoal;
end
L=sceneData.L;

%default values
x0=[L/2;L/2];
TFinal=50;
k=1;
modelName='direct';
controlArgs={};
deltaDisturbance=0;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'x0'
            ivarargin=ivarargin+1;
            x0=varargin{ivarargin};
        case 'tfinal'
            ivarargin=ivarargin+1;
            TFinal=varargin{ivarargin};
        case 'model'
            ivarargin=ivarargin+1;
            modelName=varargin{ivarargin};
        case 'controlargs'
            ivarargin=ivarargin+1;
            controlArgs=[controlArgs varargin{ivarargin}];
        case 'disturbance'
            ivarargin=ivarargin+1;
            disturbance=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

%setup model to use
switch lower(modelName)
    case 'direct'
        model=@modelDirect;
        controlArgs=[1 controlArgs];
    case 'unicycle'
        model=@modelCart;
        controlArgs=[controlArgs {'dmax',2}];
    case 'integral'
        model=@modelIntegral;
end

%compute desired bearings at goal location
YGoal=@(t) bearingCompute(XGoal(t),XLandmarks);
%derivatie term testing
%XGoalDot=@(t) functionalDerivative(XGoal(t),t);
XGoalDot=@(t) [-0.1*sin(0.1*t);0.1*cos(0.1*t)];
YGoalDot=@(t) bearingComputeDerivative(XGoalDot(t),YGoal(t),ones(1,size(YGoal(t),2))); %last input ny: ranges corresponding to the bearings?
switch lower(modelName)
    case 'integral'
        %z0=[x0; zeros(2,1)];
        %[t,x]=ode45(@(t,x) model(x,XLandmarksEval,YGoal(t),YGoalDot(t),XGoalDot(t),funs,controlArgs), [0 TFinal], z0);
        x=x';
    otherwise
        %use ode45 to integrate the closed-loop ODE
        [t,x]=ode45(@(t,x) model(x,XLandmarksEval,YGoal(t),funs,controlArgs), [0 TFinal], x0);
        x=x';
end

%closed-loop ODE for single integrator
function dx=modelDirect(XEval,X,YGoal,funs,controlArgs)
if ~iscell(controlArgs)
    controlArgs={controlArgs};
end
%compute measurements
YEval=bearingCompute(XEval,X);
%evaluate control law
dx=controlArgs{1}*bearingControlDirect(YEval,YGoal,funs,controlArgs{2:end});

%closed-loop ODE for kinematic cart (unicycle)
function dx=modelCart(XEval,X,YGoal,funs,controlArgs)
%compute measurements
YEval=bearingCompute(XEval(1:2),X);
%evaluate control law
u=cartBearingControl(XEval(3),YEval,YGoal,funs,controlArgs{:});
%apply control law to model
dx=cartModel(XEval,u);

%closed-loop ODE for integral-direct control
function dz=modelIntegral(ZEval,X,YGoal,YGoalDot,XGoalDot,funs,controlArgs)
%compute measurements
s=size(ZEval,1)/2;
YEval=bearingCompute(ZEval(1:s),X);
%evaluate control law
dz=bearingControlDirectIntegral(ZEval,YEval,YGoal,funs);
%apply derivative term
%dy=bearingCompute(ZEval(1:s)+dz(1:s),X)-YEval;
%xDot=bearingComputeDerivativeFeedback(dy,YEval,ones(1,size(YEval,2)));
%dz(1:s,:)=dz(1:s,:)+XGoalDot+0.5*xDot;
XDot=bearingComputeDerivativeFeedback(YGoalDot,YEval,ones(1,size(YGoal,2)));
dz(1:s,:)=dz(1:s,:)+1*XDot;
