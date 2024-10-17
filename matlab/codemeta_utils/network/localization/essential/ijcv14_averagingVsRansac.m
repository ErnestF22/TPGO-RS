function errors=ijcv14_averagingVsRansac(data,iMatch)
NSamples=2000;
NSampleAverage=50;
threshold=3e-3;
methodResiduals='sampsonabs';
fileNameLog=[mfilename '_log'];

[x1,x2]=sfm_getFeatureLocationFromMatchId(data,iMatch,'normalized','member','matchFiltered');
x=cat(3,x1,x2);

disp('### RANSAC Sampling')
iSampling=1;

samplingText{iSampling}='Best Of Batch';
[E{iSampling}, residuals{iSampling}, flagInlier{iSampling},allResiduals{iSampling},ESamples{iSampling},ESamplesRansac{iSampling}] = ...
    sfm_rawEssentialRansac(x1,x2,NSamples,threshold,'progressbar','keepBestSample');
iSampling=iSampling+1;

samplingText{iSampling}='Validated';
[E{iSampling}, residuals{iSampling}, flagInlier{iSampling},allResiduals{iSampling},ESamples{iSampling},ESamplesRansac{iSampling}] = ...
    sfm_rawEssentialRansac(x1,x2,NSamples,threshold,'progressbar','methodE','5pt+validation');
iSampling=iSampling+1;

%Note: residuals, flagInlier and allResiduals are just stored in the log,
%and never actually used otherwise
save(fileNameLog)

%RAverage and TAverage are cell arrays where each cell contains the
%rotations and (optionally) translations produced by each method after each
%RANSAC samples. legendTextR, legendTextT contain a descriptive string for
%each method. Note that some methods do not produce translations, so
%length(RAverage) is not equal to length(TAverage)
data.RAverage={};
data.TAverage={};
data.QAverage={};
data.legendTextR={};
data.legendTextT={};
data.legendTextQ={};

[RTruth,TTruth]=G2RT(G2GNormalized(data.matchFilteredPoseTruth(:,:,iMatch)));
QTruth=essential_fromRT(RTruth,TTruth);

NSamplings=iSampling-1;
for iSampling=1:NSamplings
    disp('### Solve twisted pair ambiguity')
    NPointsTwistedPair=min(20,size(x1,2));
    [RSamples,TSamples]=epipolarEToRT(ESamples{iSampling},x1(:,1:NPointsTwistedPair),x2(:,1:NPointsTwistedPair),'progressbar');
    [RRansac,TRansac]=epipolarEToRT(E{iSampling},x1(:,1:NPointsTwistedPair),x2(:,1:NPointsTwistedPair));
    [RSamplesRansac,TSamplesRansac]=epipolarEToRT(ESamplesRansac{iSampling}(:,:,1:NSampleAverage),x1(:,1:NPointsTwistedPair),x2(:,1:NPointsTwistedPair));
    QSamples=essential_fromRT(RSamples,TSamples);
    QRansac=essential_fromRT(RRansac,TRansac);
    QSamplesRansac=essential_fromRT(RSamplesRansac,TSamplesRansac);

    disp('### Pose extraction methods')
    testText='Rotation averaging L2';
    RAverageR=cumAverageSamplesR(RSamples,NSampleAverage,2);
    data=accumulateData(data,testText,samplingText{iSampling},RAverageR);

    testText='Rotation averaging L1';
    RAverageR=cumAverageSamplesR(RSamples,NSampleAverage,1);
    data=accumulateData(data,testText,samplingText{iSampling},RAverageR);

    testText='Essential averaging L2';
    [RAverageQ,TAverageQ,QAverageQ]=cumAverageSamplesQ(QSamples,NSampleAverage,2);
    data=accumulateData(data,testText,samplingText{iSampling},RAverageQ,TAverageQ,QAverageQ);
    
    testText='Essential averaging L1';
    [RAverageQ,TAverageQ,QAverageQ]=cumAverageSamplesQ(QSamples,NSampleAverage,1);
    data=accumulateData(data,testText,samplingText{iSampling},RAverageQ,TAverageQ,QAverageQ);

    testText='Essential averaging L1 - Sampson';
    [RRefined,TRefined,QRefined]=essentialRefine(QAverageQ,x,threshold,methodResiduals);
    data=accumulateData(data,testText,samplingText{iSampling},RRefined,TRefined,QRefined);

    testText='Samples';
    NSampleAverageActual=min(NSampleAverage,size(RSamples,3));
    data=accumulateData(data,testText,samplingText{iSampling},...
        RSamples(:,:,1:NSampleAverageActual),TSamples(:,1:NSampleAverageActual),QSamples(:,:,1:NSampleAverageActual));

    testText='Samples RANSAC';
    data=accumulateData(data,testText,samplingText{iSampling},RSamplesRansac,TSamplesRansac,QSamplesRansac);

    testText='Samples RANSAC - Sampson';
    [RRefined,TRefined,QRefined]=essentialRefine(QSamplesRansac,x,threshold,methodResiduals);
    data=accumulateData(data,testText,samplingText{iSampling},RRefined,TRefined,QRefined);
    
    testText='RANSAC';
    data=accumulateData(data,testText,samplingText{iSampling},...
        repmat(RRansac,[1 1 NSampleAverage]),repmat(TRansac,[1 NSampleAverage]),repmat(QRansac,[1 1 NSampleAverage]));
    
    testText='RANSAC - Sampson';
    [RRefined,TRefined,QRefined]=essentialRefine(QRansac,x,threshold,methodResiduals);
    data=accumulateData(data,testText,samplingText{iSampling},...
        repmat(RRefined,[1 1 NSampleAverage]),repmat(TRefined,[1 NSampleAverage]),repmat(QRefined,[1 1 NSampleAverage]));

    save(fileNameLog)
end

disp('### Compute error measures')
NMethodsR=length(data.RAverage);
errR=zeros(NSampleAverage,NMethodsR);
for iMethod=1:NMethodsR
    errR(:,iMethod)=padrow(rot_dist(RTruth,data.RAverage{iMethod}),NSampleAverage);
end

NMethodsT=length(data.TAverage);
errT=zeros(NSampleAverage,NMethodsT);
for iMethod=1:NMethodsT
    errT(:,iMethod)=padrow(sphere_dist(TTruth,data.TAverage{iMethod}),NSampleAverage);
end

NMethodsQ=length(data.QAverage);
errQ=zeros(NSampleAverage,NMethodsQ);
for iMethod=1:NMethodsQ
    errQ(:,iMethod)=padrow(essential_dist(QTruth,data.QAverage{iMethod},'flagUseMex',false),NSampleAverage);
end

save(fileNameLog)

errors.errR=errR;
errors.errT=errT;
errors.errQ=errQ;
errors.legendTextR=data.legendTextR;
errors.legendTextT=data.legendTextT;
errors.legendTextQ=data.legendTextQ;

function r=padrow(r,l)
r=shiftdim(r)';
lr=size(r,2);
if lr<l
    r=[r NaN(1,l-lr)];
end

%Cumulative average of the sampled rotations
function RAverage=cumAverageSamplesR(RSamples,NSamples,LNorm)
RAverage=cumAverageSamples(rot_funs,RSamples,NSamples,LNorm);

%Cumulative average of the sampled essential matrices on the essential
%manifold
function [RAverage,TAverage,QAverage]=cumAverageSamplesQ(QSamples,NSamples,LNorm)
lf=essential_signed_funs();
lf.log=@(Q1,Q2,varargin) essential_log(Q1,Q2,'signed','flagusemex',false,varargin{:});
lf.dist=@(Q1,Q2,varargin) essential_dist(Q1,Q2,'signed','flagusemex',false,varargin{:});
QAverage=cumAverageSamples(lf,QSamples,NSamples,LNorm);
[RAverage,TAverage]=essential_toRT(QAverage);

%Cumulative average of points using the Weiszfeld algorithm
function yAverage=cumAverageSamples(lf,ySamples,NSamples,LNorm)
NSamples=min(NSamples,size(ySamples,3));
flagDisplaIt=false;
yAverage=zeros(size(ySamples,1),size(ySamples,2),NSamples);
w=getTextWaitBar(NSamples);
optsLieMinimize={'norm',LNorm,'maxit',30};
for iSample=1:NSamples
    yInit=lie_minimize_tangentAverage_getInitialPoint(lf,ySamples(:,:,1:iSample),'method','twominaverage');
    if iSample==NSamples && flagDisplaIt
        optsLieMinimize=[optsLieMinimize 'displayIt'];
    end
    yAverage(:,:,iSample)=lie_minimize_tangentAverage(lf,yInit,ySamples(:,:,1:iSample),...
        optsLieMinimize{:});
    w(iSample)
end

%starting from Q, find inliers and then refine Sampson error
function [RRefined,TRefined,QRefined]=essentialRefine(Q,x,threshold,methodResiduals)
E=essential_toE(Q);
residuals=epipolarResiduals(E,x(:,:,1),x(:,:,2),methodResiduals);
inliers=residuals<threshold;
QRefined=zeros(size(Q));
for iQ=1:size(Q,3)
    QRefined(:,:,iQ)=essential_minimizeEpipolarCostSampsonSq(...
        x(:,inliers(iQ,:),:),Q(:,:,iQ));
end
[RRefined,TRefined]=essential_toRT(QRefined);

%function to accumulate results from different methods
function data=accumulateData(data,testText,samplingText,R,T,Q)
testText=[testText ' - ' samplingText];
disp(['^--' testText])
data.RAverage=[data.RAverage R];
data.legendTextR=[data.legendTextR testText];
if exist('T','var')
    data.TAverage=[data.TAverage T];
    data.legendTextT=[data.legendTextT testText];
end
if exist('Q','var')
    data.QAverage=[data.QAverage Q];
    data.legendTextQ=[data.legendTextQ testText];
end
