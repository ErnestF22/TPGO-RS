function data=POCsfm_matchExtractQuickMatch(data,optsQuickMatch)
matchMemberName='match_quickMatch';
NImages=sfm_getImageNumber(data);
ratioDensity=0.25;
ratioInterCluster=0.75;
%ratioInterCluster=1;
ratioIntraCluster=1;
thresholdBreakTree=Inf;

param = [{'ratioDensity',ratioDensity,...
    'ratioInterCluster',ratioInterCluster,...
    'ratioIntraCluster',ratioIntraCluster,...
    'threshold',thresholdBreakTree} optsQuickMatch];

[dataMatch]=sfm_featureExtractJointMatchData(data,'NImages',NImages);
display('Distance matrix computation')
D=sqrt(euclideanDistMatrix(double(dataMatch.descriptor)));
display('QuickMatch')
[membershipMatches,dataMatch.QuickMatchInfo]=quickshift_matching(D,dataMatch.membershipPrior,'gaussian','densityLogAmplify',param{:});
membershipPrior=dataMatch.membershipPrior;

%Multimatch clusters in the form of a matrix of indicators
X = sqrt(euclideanDistMatrix(membershipMatches));
X = X == 0;

%Store information from multimatch clusters in X into the data.match structure
NMatchPairs=NImages*(NImages-1)/2;
match=repmat(struct('idxImg',[],'idxMatch',[]),1,NMatchPairs);
cntMatch=1;
for iImg = 1:NImages-1
    for jImg = iImg+1:NImages
        match(cntMatch).idxImg=[iImg;jImg];
        flagInd1 = membershipPrior==iImg;
        flagInd2 = membershipPrior==jImg;
        Y = X(flagInd1,flagInd2);
        [id1,id2] = find(Y);
        match(cntMatch).idxMatch = [id1,id2]';
        cntMatch=cntMatch+1;
    end
end
data.(matchMemberName)=match;