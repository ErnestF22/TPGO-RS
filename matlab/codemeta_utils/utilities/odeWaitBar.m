function status=odeWaitBar(t,varargin)
persistent h t0

status=0;
if isempty(h) || length(t)>1
    t0=t(1);
    h=getTextWaitBar(t(2)-t0);
    h(0);
elseif isempty(t)
    fprintf('\n')
else
    h(t-t0)
end
