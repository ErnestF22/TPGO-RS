function POCAveragingVsRansac
%resetRands()
fileNameCache=[mfilename '_cache.mat'];
NSamples=100;
NSampleAverage=40;
threshold=3e-3;
iMatch=11;

nFig=1;
flagDisplayRotErrDist=false;

load('cvpr13_fountain_data','data')
[x1,x2]=sfm_getFeatureLocationFromMatchId(data,iMatch,'normalized');

% if exist(fileNameCache,'file')
%     load(fileNameCache)
%     disp('Loaded cache')
% else

disp('RANSAC Sampling')
[E, ~, ~,~,ESamples] = ...
    sfm_rawEssentialRansac8pt(x1,x2,NSamples,threshold,'progressbar');

disp('Twisted pair ambiguity')
NPointsTwistedPair=10;
[RSamples,TSamples]=epipolarEToRT(ESamples,x1(:,1:NPointsTwistedPair),x2(:,1:NPointsTwistedPair),'progressbar');

save(fileNameCache,'E','ESamples','RSamples','TSamples')
% disp('Saved cache')
% end

[RRansac,TRansac]=epipolarEToRT(E,x1(:,1:NPointsTwistedPair),x2(:,1:NPointsTwistedPair));

[RTruth,TTruth]=G2RT(G2GNormalized(data.matchPoseTruth(:,:,iMatch)));

if flagDisplayRotErrDist
    figure(nFig); nFig=nFig+1;
    cumDistPerc(rot_dist(RTruth,RSamples));
    hold on
    plot(rot_dist(RTruth,RRansac)*[1 1],[0 100],'r:')
    hold off
    title('Rotation errors for each trial')
end

RAverage={};
TAverage={};
legendTextR={};
legendTextT={};

testText='Rotation averaging L2';
disp(testText)
RAverageR=cumAverageSamplesR(RSamples,NSampleAverage,2);
RAverage=[RAverage RAverageR];
legendTextR=[legendTextR testText];

testText='Rotation averaging L1';
disp(testText)
RAverageR=cumAverageSamplesR(RSamples,NSampleAverage,1);
RAverage=[RAverage RAverageR];
legendTextR=[legendTextR testText];

testText='Essential averaging L2';
disp(testText)
QSamples=essential_fromRT(RSamples,TSamples);
[RAverageQ,TAverageQ]=cumAverageSamplesQ(QSamples,NSampleAverage,2);
RAverage=[RAverage RAverageQ];
TAverage=[TAverage TAverageQ];
legendTextR=[legendTextR testText];
legendTextT=[legendTextT testText];

testText='Essential averaging L1';
disp(testText)
QSamples=essential_fromRT(RSamples,TSamples);
[RAverageQ,TAverageQ]=cumAverageSamplesQ(QSamples,NSampleAverage,1);
RAverage=[RAverage RAverageQ];
TAverage=[TAverage TAverageQ];
legendTextR=[legendTextR testText];
legendTextT=[legendTextT testText];

NMethodsR=length(RAverage);
errR=zeros(NSampleAverage,NMethodsR);
for iMethod=1:NMethodsR
    errR(:,iMethod)=rot_dist(RTruth,RAverage{iMethod});
end

NMethodsT=length(TAverage);
errT=zeros(NSampleAverage,NMethodsT);
for iMethod=1:NMethodsT
    errT(:,iMethod)=sphere_dist(TTruth,TAverage{iMethod});
end
    
figure(nFig); nFig=nFig+1;
plot(errR)
hold on
plot(rot_dist(RTruth,RSamples(:,:,1:NSampleAverage)),'r-.')
plot([1 NSampleAverage],rot_dist(RTruth,RRansac)*[1 1],'r:')
hold off
legend(legendTextR{:},'Samples','RANSAC solution')
title('Rotation errors')

figure(nFig); nFig=nFig+1;
plot(errT)
hold on
plot(sphere_dist(TTruth,TSamples(:,1:NSampleAverage)),'r-.')
plot([1 NSampleAverage],sphere_dist(TTruth,TRansac)*[1 1],'r:')
hold off
legend(legendTextT{:},'Samples','RANSAC solution')
title('Translation angle')


function plotRotErrorsAverage(RTruth,RAverage,RSamples,RRansac)
NSampleAverage=size(RAverage,3);

plot(rot_dist(RTruth,RAverage))
hold on
plot(rot_dist(RTruth,RSamples(:,:,1:NSampleAverage)),'r')
plot([1 NSampleAverage],rot_dist(RTruth,RRansac)*[1 1],'r:')
hold off
legend('Errors with averaging','Errors of samples','Errors of RANSAC solution')


function RAverage=cumAverageSamplesR(RSamples,NSamples,LNorm)
RAverage=cumAverageSamples(rot_funs,RSamples,NSamples,LNorm);

function [RAverage,TAverage]=cumAverageSamplesQ(QSamples,NSamples,LNorm)
QAverage=cumAverageSamples(essential_signed_funs,QSamples,NSamples,LNorm);
[RAverage,TAverage]=essential_toRT(QAverage);

function yAverage=cumAverageSamples(lf,ySamples,NSamples,LNorm)
yAverage=zeros(size(ySamples,1),size(ySamples,2),NSamples);
w=getTextWaitBar(NSamples);
for iSample=1:NSamples
    yInit=lie_minimize_tangentAverage_getInitialPoint(lf,ySamples(:,:,1:iSample),'method','twominaverage');
    yAverage(:,:,iSample)=lie_minimize_tangentAverage(lf,yInit,ySamples(:,:,1:iSample),...
        'norm',LNorm,'maxit',30);
    w(iSample)
end

