function POCPlotSphereMap

A=randn(3);
f=@(x) A*x;

NRadii=4;
NLatitudes=5;
NLongitudes=8;
stepRadii=1;

p=zeros(NRadii,NLatitudes,NLongitudes,3);
for iRadii=1:NRadii
    r=stepRadii*(iRadii-1);
    for iLongitudes=1:NLongitudes
        theta=2*pi*(iLongitudes-1)/NLongitudes;
        for iLatitudes=1:NLatitudes
            phi=pi*(iLatitudes-1)/(NLatitudes-1);
            p(iRadii,iLatitudes,iLongitudes,1)=r*sin(phi)*cos(theta);
            p(iRadii,iLatitudes,iLongitudes,2)=r*sin(phi)*sin(theta);
            p(iRadii,iLatitudes,iLongitudes,3)=r*cos(phi);
        end
    end
end

p=evalfunVec(f,p);

cmap=winter(NRadii);
for iRadii=1:NRadii
    plotArgs={'b-','color',cmap(iRadii,:)};
    c=cmap(iRadii,:);
    if iRadii~=NRadii
        pPlot=reshape(p(iRadii:iRadii+1,:,:,:),2,[],3);
        plot3(pPlot(:,:,1),pPlot(:,:,2),pPlot(:,:,3),plotArgs{:})
    end
    hold on
    for iLatitudes=1:NLatitudes-1
        pPlot=reshape(permute(p(iRadii,iLatitudes:iLatitudes+1,:,:),[2 1 3 4]),2,[],3);
        plot3(pPlot(:,:,1),pPlot(:,:,2),pPlot(:,:,3),plotArgs{:})
    end
    for iLongitudes=1:NLongitudes
        idxLongitudes=mod((iLongitudes:iLongitudes+1)-1,NLongitudes)+1;
        pPlot=reshape(permute(p(iRadii,:,idxLongitudes,:),[3 1 2 4]),2,[],3);
        plot3(pPlot(:,:,1),pPlot(:,:,2),pPlot(:,:,3),plotArgs{:})
    end    
end
hold off
axis equal
