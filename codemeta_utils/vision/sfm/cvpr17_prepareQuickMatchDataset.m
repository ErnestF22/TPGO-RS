function cvpr17_prepareQuickMatchDataset
datasetName='fountain_dense';
fileNameData=[mfilename '_data.mat'];
NImages=3;
if ~exist(fileNameData,'file')
    data=sfm_initDataFromDir(['/Users/tron/Documents/UPenn/datasets/' datasetName]);
    data=sfm_featureExtract(data,'flagAutoFirstOctave',false,'optsSift',{'FirstOctave', 2},'showStats');
    %data=sfm_featureExtract(data,'showStats');
    data=sfm_addCalibration(data);
    data=sfm_featureNormalize(data);
    G=cvlabLoadCameras(['/Users/tron/Documents/UPenn/datasets/' datasetName]);
    data=sfm_poseAdd(data,G);
    save(fileNameData,'data')
else
    load(fileNameData)
end
pMatch=runQuickMatchSfM(data,NImages);
data.NImages=NImages;
%pMatch=runQuickMatchSfM(data);
pMatch=addEssentialToMatches(data,pMatch);
pMatchPair=pairwiseMatches(data);
pMatchPair=addEssentialToMatches(data,pMatchPair);
thresholds=linspace(0,0.5,30);
meanacc.quickshift=computeMeanAccuracy(pMatch,thresholds);
meanacc.pair=computeMeanAccuracy(pMatchPair,thresholds);

plot(thresholds,meanacc.quickshift,thresholds,meanacc.pair)
legend('quickshift','pair')
%keyboard

function [meanacc,pMatch]=computeMeanAccuracy(pMatch,thresholds)
NImages=size(pMatch,2);
for iImage=1:NImages
    for jImage=(iImage+1):NImages
        %displayMatch(data,pMatchPair,iImage,jImage)
        %disp(pMatch(iImage,jImage).Ematrix)
        r=evaluateMatch(pMatch,iImage,jImage);
        pMatch(iImage,jImage).acc=computeAccuracy(r,thresholds,max(pMatch(iImage,jImage).nFeature));
    end
end
meanacc=mean(cat(1,pMatch.acc));

function acc=computeAccuracy(r,thresholds,NTotal)
NThresholds=length(thresholds);
acc=zeros(1,NThresholds);
for iThreshold=1:NThresholds
    acc(iThreshold)=sum(r<thresholds(iThreshold))/NTotal;
end

function pMatch=pairwiseMatches(data)
disp('Pairwise matches')
tic
NImages=data.NImages;
for iImage=1:NImages
    for jImage=(iImage+1):NImages
        feature1=data.feature(iImage);
        feature2=data.feature(jImage);
        desc1=feature1.descriptor;
        desc2=feature2.descriptor;
        m12=vl_ubcmatch(desc1,desc2);
        idx1=m12(1,:);
        idx2=m12(2,:);
        pMatch(iImage,jImage).location={feature1.location(:,idx1) feature2.location(:,idx2)};
        pMatch(iImage,jImage).locationNormalized={feature1.locationNormalized(:,idx1) feature2.locationNormalized(:,idx2)};
        pMatch(iImage,jImage).nFeature=[size(desc1,2) size(desc2,2)];
    end
end
toc

function pMatch=addEssentialToMatches(data,pMatch)
NImages=size(pMatch,2);
for iImage=1:NImages
    for jImage=(iImage+1):NImages
        Gi=data.poseTruth(:,:,iImage);
        Gj=data.poseTruth(:,:,jImage);
        Eij=epipolarBuildEFromG(Gi,Gj,'references');
        Ki=data.calibration(:,:,iImage);
        Kj=data.calibration(:,:,jImage);
        pMatch(iImage,jImage).Ematrix=Eij;
    end
end

function r=evaluateMatch(pMatch,iImage,jImage)
E=pMatch(iImage,jImage).Ematrix;
x1=pMatch(iImage,jImage).locationNormalized{1};
x2=pMatch(iImage,jImage).locationNormalized{2};
r=epipolarResiduals(E',x1,x2,'linedistanceabs');

function displayMatch(data,pMatch,iImage,jImage)
imgi=sfm_getImageById(data,iImage);
imgj=sfm_getImageById(data,jImage);
xi=pMatch(iImage,jImage).location{1};
xj=pMatch(iImage,jImage).location{2};
sfm_rawDisplayMatch(imgi,imgj,xi,xj)
