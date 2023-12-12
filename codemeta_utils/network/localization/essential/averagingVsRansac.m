function errors=averagingVsRansac(data,iMatch)
NSamples=2000;
NSampleAverage=50;
threshold=3e-3;

[x1,x2]=sfm_getFeatureLocationFromMatchId(data,iMatch,'normalized');

disp('RANSAC Sampling')
[E, ~, ~,~,ESamples,ESamplesRansac] = ...
    sfm_rawEssentialRansac(x1,x2,NSamples,threshold,'progressbar','keepBestSample');

disp('Solve twisted pair ambiguity')
NPointsTwistedPair=min(20,size(x1,2));
[RSamples,TSamples]=epipolarEToRT(ESamples,x1(:,1:NPointsTwistedPair),x2(:,1:NPointsTwistedPair),'progressbar');
[RRansac,TRansac]=epipolarEToRT(E,x1(:,1:NPointsTwistedPair),x2(:,1:NPointsTwistedPair));
[RSamplesRansac,TSamplesRansac]=epipolarEToRT(ESamplesRansac(:,:,1:NSampleAverage),x1(:,1:NPointsTwistedPair),x2(:,1:NPointsTwistedPair));

[RTruth,TTruth]=G2RT(G2GNormalized(data.matchPoseTruth(:,:,iMatch)));

%RAverage and TAverage are cell arrays where each cell contains the
%rotations and (optionally) translations produced by each method after each
%RANSAC samples. legendTextR, legendTextT contain a descriptive string for
%each method. Note that some methods do not produce translations, so
%length(RAverage) is not equal to length(TAverage)
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

testText='Samples';
disp(testText)
RAverage=[RAverage RSamples(:,:,1:NSampleAverage)];
TAverage=[TAverage TSamples(:,1:NSampleAverage)];
legendTextR=[legendTextR testText];
legendTextT=[legendTextT testText];

testText='Samples RANSAC';
disp(testText)
RAverage=[RAverage RSamplesRansac];
TAverage=[TAverage TSamplesRansac];
legendTextR=[legendTextR testText];
legendTextT=[legendTextT testText];

testText='RANSAC';
disp(testText)
RAverage=[RAverage repmat(RRansac,[1 1 NSampleAverage])];
TAverage=[TAverage repmat(TRansac,[1 NSampleAverage])];
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

errors.errR=errR;
errors.errT=errT;
errors.legendTextR=legendTextR;
errors.legendTextT=legendTextT;

%Cumulative average of the sampled rotations
function RAverage=cumAverageSamplesR(RSamples,NSamples,LNorm)
RAverage=cumAverageSamples(rot_funs,RSamples,NSamples,LNorm);

%Cumulative average of the sampled essential matrices on the essential
%manifold
function [RAverage,TAverage]=cumAverageSamplesQ(QSamples,NSamples,LNorm)
lf=essential_signed_funs();
lf.log=@(Q1,Q2,varargin) essential_log(Q1,Q2,'signed','flagusemex',false,varargin{:});
lf.dist=@(Q1,Q2,varargin) essential_dist(Q1,Q2,'signed','flagusemex',false,varargin{:});
QAverage=cumAverageSamples(lf,QSamples,NSamples,LNorm);
[RAverage,TAverage]=essential_toRT(QAverage);

%Cumulative average of points using the Weiszfeld algorithm
function yAverage=cumAverageSamples(lf,ySamples,NSamples,LNorm)
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
