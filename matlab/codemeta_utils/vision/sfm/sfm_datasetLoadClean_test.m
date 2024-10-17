function sfm_datasetLoadClean_test
data=sfm_datasetLoadClean('fountain');
data=sfm_addMatchEssentialTruth(data);
NMatches=length(data.match);
e=NaN(1,NMatches);
for iMatch=1:NMatches
    m=data.match(iMatch).idxMatch;
    idxImg=data.match(iMatch).idxImg;
    x1=data.feature(idxImg(1)).locationNormalized(:,m(1,:));
    x2=data.feature(idxImg(2)).locationNormalized(:,m(2,:));

    e(iMatch)=max(abs(epipolarConstraintFromE(data.matchEssentialTruth(:,:,iMatch),x1,x2)));
end

disp(e)
