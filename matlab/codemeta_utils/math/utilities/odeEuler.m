%function [t,x]=odeEuler(dx,tspan,x0,options)
%Solve ODE using Euler's method (i.e., explicit, first order, fixed step
%size).
%Inputs
%   dx      Function hanlde of the form dx(t,x) describing the derivative field
%   tspan   [2 x 1] array with initial and final integration times
%   x0      Initial conditions
%   options Struct of options set with odeset (see below for which options
%           are actually used)
%
%Options used
%   'MaxStep'   Step size used for the integration
%   'OutputFnc' Function to visualize output (see odeset)
function [t,x]=odeEuler(dx,tspan,x0,options)
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


t=(tspan(1):dt:tspan(2))';
Nt=length(t);
dxTemp=dx(0,x0);
d=length(dxTemp(:));
x=zeros(Nt,d);
x(1,:)=x0;
for it=2:Nt
    dxt=dx(t(it-1),x(it-1,:)');
    x(it,:)=x(it-1,:)+(t(it)-t(it-1))*dxt(:)';
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
