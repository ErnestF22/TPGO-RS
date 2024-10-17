function iros14_plots
load('multiCalibration_testRealData_N100_data')
figDir='../../../../papers/vision/IROS14-MultiCalibration/figures';
flagSingle=false;
figDim=[260,115];
figDimBig=[260,225];
flagSaveFigure=2;

fs=setFigFontSize(8);
fn=setFigFont('Times');

figure(1)
showErrors(errorsAggregated,'cameraVelodyneRot',180/pi,flagSingle,true);
%title('Camera-Velodyne rotation')
xlabel('[deg]')
ylabel('% of total')
savefigure(fullfile(figDir,'cameraVelodyneRot'),'epsc',figDim,flagSaveFigure)

figure(2)
showErrors(errorsAggregated,'cameraVelodyneTransl',100,flagSingle);
%title('Camera-Velodyne translation')
xlabel('[cm]')
ylabel('% of total')
savefigure(fullfile(figDir,'cameraVelodyneTransl'),'epsc',figDim,flagSaveFigure)

figure(3)
showErrors(errorsAggregated,'hokuyoCamera',100,flagSingle);
%title('Hokuyo-Camera')
xlabel('[cm]')
ylabel('% of total')
savefigure(fullfile(figDir,'hokuyoCamera'),'epsc',figDim,flagSaveFigure)

figure(4)
showErrors(errorsAggregated,'hokuyoVelodyne',100,flagSingle);
%title('Hokuyo-Velodyne')
xlabel('[cm]')
ylabel('% of total')
savefigure(fullfile(figDir,'hokuyoVelodyne'),'epsc',figDim,flagSaveFigure)


figure(5)
GComp=RT2G(rot([0;0;pi])*rot([-pi/2;0;0]),zeros(3,1))*invg(poses{end}.GCamera);
methodAbsolutePoses='reference';
draw3dcameraFromG(GComp*poses{end}.GCamera,'scale',0.2,'methodAbsolutePoses',methodAbsolutePoses,'color2',[0.1 0.8 0.1],'color1',[0.1 0.8 0.1]);
hold on
draw3dcameraFromG(GComp*poses{end}.GVelodyne,'shape','velodyne','scale',0.2,'methodAbsolutePoses',methodAbsolutePoses,'alpha',1,'color1',[0.8 0.1 0.1]);
draw3dcameraFromG(GComp*poses{end}.GHokuyo,'shape','hokuyo','scale',0.25,'methodAbsolutePoses',methodAbsolutePoses,'color1',[0.7 0.7 0.9]);
hold off
view(-56,30)
set(gcf,'Renderer','Painters')
axis equal
savefigure(fullfile(figDir,'poses'),'epsc',figDimBig,flagSaveFigure)

setFigFontSize(fs);
setFigFont(fn);


function showErrors(errorsAggregated,fieldName,factor,flagSingle,flagLegend)
if ~exist('flagLegend','var')
    flagLegend=false;
end
if flagSingle
    e=errorsAggregated.(fieldName);
else
    e=[ errorsAggregated.([fieldName 'Measured'])...
        errorsAggregated.(fieldName)...
        errorsAggregated.([fieldName 'Indirect'])...
        errorsAggregated.([fieldName 'Hand'])];
end
e=abs(e)*factor;
cumDistBoxPerc(e)
if ~flagSingle && flagLegend
    legend('Initial','Averaged','Indirect','Manual')
end
disp(['# ' fieldName ' median'])
disp(median(e,1))
