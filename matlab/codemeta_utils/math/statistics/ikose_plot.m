function ikose_plot(r,varargin)
s=ikose(r,varargin{:});
cumDistPerc(r)
hold on
plot([s s],[0 100],'k')
hold off
