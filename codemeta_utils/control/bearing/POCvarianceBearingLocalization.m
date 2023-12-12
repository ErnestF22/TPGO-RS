function POCvarianceBearingLocalization
resetRands();
nLandmarkSet=1;
fileNameSave=[mfilename '_' num2str(nLandmarkSet) '_data'];
funs=bearingCostFunctions('cosine');
NX=30;
switch nLandmarkSet
    case 1
        XLandmarks0=2*(rand(2,NX)-0.5);
    case 2
        XLandmarks0=[-1 -1 1 1; -1 1 -1 1];
    case 3
        XLandmarks0=[0 sin(2/3*pi) sin(4/3*pi); 1 cos(2/3*pi) cos(4/3*pi)];
end

ND=5;
z=cell(1,ND);
for iD=1:ND
    D=diag([1 1/iD]);
    XLandmarks=D*XLandmarks0;
    %XLandmarks=2*D*[linspace(-0.5,0.5,NX); rand(1,NX)-0.5];
    xmax=3;%2*max(diag(D));
    x=linspace(-xmax,xmax(1),30);
    vContour=linspace(0,10000,100);

    [xx,yy]=meshgrid(x,x);
    z{iD}=evalfunVec(@(x) variance(x,XLandmarks,funs),cat(3,xx,yy));
    figure(iD)
    contour(xx,yy,z{iD},vContour);
    disp(min(z{iD}(:)))
    disp(max(z{iD}(:)))
    hold on
    plot(XLandmarks(1,:),XLandmarks(2,:),'r*');
    hold off
    axis equal
    axis tight
end
save(fileNameSave)


function v=variance(XEval,XLandmarks,funs)
[YEval,nYEval]=bearingCompute(XEval,XLandmarks);
v=trace(bearingCostGeneral_covarianceEstimate(YEval,YEval,nYEval,funs));
