%Integrate second order ODE on arbitrary spaces using Euler's method
%function [t,x]=generic_odeEuler(dx,tspan,x0,funTransition,options)
%Inputs
%   funTransition   Handle to a function with prototype
%       xNew=funTransition(x,dt,dx)
function [t,x]=generic_odeEuler(dx,tspan,x0,funTransition,options)
if ~exist('options','var')
    options=[];
end
dt=odeget(options,'MaxStep',0.01);

% Handle the output function
if nargout > 0
  outputFcn = odeget(options,'OutputFcn',[]);
else
  outputFcn = odeget(options,'OutputFcn',@odeplot);
end

if isempty(outputFcn)
  flagOutputFcn = false;
else
  flagOutputFcn = true;
  outputIdx = odeget(options,'OutputSel',1:length(x0),'fast');
end

% Initialize the output function.
if flagOutputFcn
  feval(outputFcn,tspan,x0(outputIdx),'init');
end

%Initialization
t=(tspan(1):dt:tspan(2))';
Nt=length(t);
dxTemp=dx(0,x0);
d=length(dxTemp(:));
x=zeros(Nt,d);
x(1,:)=x0;
%[R,w]=rotDyn_stateUnpack(x0);

%Step by step integration
for it=2:Nt
    xCurrent=x(it-1,:)';
    dxCurrent=dx(t(it-1),xCurrent);
    dtCurrent=(t(it)-t(it-1));
    x(it,:)=funTransition(xCurrent,dtCurrent,dxCurrent);

    %Call output function if provided
    if flagOutputFcn
      stop=feval(outputFcn,t(it),x(it,outputIdx)','');
      if stop
          break
      end
    end
end

% Cleanup for the output function.
if flagOutputFcn
  feval(outputFcn,[],[],'done');
end

