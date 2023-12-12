function icra13_varianceBearingLocalization
resetRands();
nLandmarkSet=2;
fileNameSave=[mfilename '_' num2str(nLandmarkSet) '_data'];
funs=bearingCostFunctions('cosine');
NX=30;
NGrid=200;
switch nLandmarkSet
    case 1
        XLandmarks0=6*(rand(2,NX)-0.5);
    case 2
        XLandmarks0=3*[-1 -1 1 1; -1 1 -1 1];
    case 3
        XLandmarks0=3*[0 sin(2/3*pi) sin(4/3*pi); 1 cos(2/3*pi) cos(4/3*pi)];
    case 4
        load('icra13_singleTrajectoryComparison_data','sceneData');
        XLandmarks0=sceneData.XLandmarks;
end

ND=5;
z=cell(1,ND);
u=cell(1,ND);
v=cell(1,ND);
XLandmarks=cell(1,ND);
for iD=1:ND
    disp(['Set ' num2str(iD) '/' num2str(ND)])
    D=diag([1 1/iD]);
    XLandmarks{iD}=D*XLandmarks0;
    %XLandmarks=2*D*[linspace(-0.5,0.5,NX); rand(1,NX)-0.5];
    xmax=10;%2*max(diag(D));
    x=linspace(-xmax,xmax(1),NGrid);
    vContour=linspace(0,10000,100);

    [xx,yy]=meshgrid(x,x);
    z{iD}=zeros(size(xx));
    u{iD}=zeros([size(xx) 2]);
    v{iD}=zeros([size(xx) 2]);
    w=getTextWaitBar(NGrid);
    w(0)
    for ix=1:NGrid
        for iy=1:NGrid
            [z{iD}(ix,iy),u{iD}(ix,iy,:),v{iD}(ix,iy,:)]=variance([xx(ix,iy);yy(ix,iy)],XLandmarks{iD},funs);
        end
        w(ix)
    end
    %z{iD}=evalfunVec(@(x) variance(x,XLandmarks{iD},funs),cat(3,xx,yy));
    figure(iD)
    contour(xx,yy,z{iD},vContour);
    disp(min(z{iD}(:)))
    disp(max(z{iD}(:)))
    hold on
    plot(XLandmarks{iD}(1,:),XLandmarks{iD}(2,:),'r*');
    hold off
    axis equal
    axis tight
end
save(fileNameSave)


function [sigma,u,v]=variance(XEval,XLandmarks,funs)
[YEval,nYEval]=bearingCompute(XEval,XLandmarks);
Sigma=bearingCostGeneral_covarianceEstimate(YEval,YEval,nYEval,funs);
sigma=trace(Sigma);
[U,S,V]=svd(Sigma);
US=U*S;
u=U(:,1);
v=U(:,2);

