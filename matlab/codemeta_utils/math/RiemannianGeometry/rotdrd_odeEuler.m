function [t,x]=rotdrd_odeEuler(dx,tspan,x0,options)
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
[R,w]=rotDyn_stateUnpack(x0);
for it=2:Nt
    dxCurrent=dx(t(it-1),x(it-1,:)');
    dtCurrent=(t(it)-t(it-1));
    [~,dw]=rotDyn_stateUnpack(dxCurrent);
    R=R*rot(dtCurrent*w);
    w=w+dtCurrent*dw;
    x(it,:)=rotDyn_statePack(R,w);
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