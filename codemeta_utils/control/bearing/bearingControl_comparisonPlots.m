function bearingControl_comparisonPlots

funsBearings=bearingCostFunctions('cosine');
funsAngle.f=@(x) x^2/2;

Nr=5;
Ntheta=30;
x=reshape(permute(polarGrid(Nr,Ntheta),[3 1 2]),2,[]);
Nx=size(x,2);
xLandmark=[0;0];

y0=bearingCompute([1;0],xLandmark);
y=bearingCompute(x,xLandmark);
dxGradient=bearingControlDirect(y,y0,funsBearings);
dxGeodesic=bearingField(y,y0,'geodesic',funsAngle);
dxProjection=bearingField(y,y0,'projection');

plotField(x,dxGradient,'b')
hold on
plotField(x,dxGeodesic,'r')
plotField(x,dxProjection,'c')
hold off



