function setMyColorOrder(N)
cmap=[
         0         0    1.0000
         0    0.5000         0
    1.0000         0         0
         0    0.7500    0.7500
    0.7500         0    0.7500
    0.7500    0.7500         0
    0.2500    0.2500    0.2500
]*0.8+0.2;
if ~exist('N','var')
    N=size(cmap,1);
end

cmap=cmap(1:N,:);

set(gca,'ColorOrder',cmap,'LineStyleOrder',{'-','--',':'},'NextPlot','replaceChildren');
