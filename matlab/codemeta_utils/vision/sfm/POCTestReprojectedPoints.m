function POCTestReprojectedPoints
load('sfmdata_fountain')
m=data.match(1).idxMatch;
idxImg=data.match(1).idxImg;

structId=zeros(size(m));
for k=1:2
    structId(k,:)=data.feature(idxImg(k)).structureFilteredMembership(m(k,:));
end
flagM=all(structId>0) & structId(1,:)==structId(2,:);
m=m(:,flagM);

x1=data.feature(idxImg(1)).locationReprojectedNormalized(:,m(1,:));
x2=data.feature(idxImg(2)).locationReprojectedNormalized(:,m(2,:));

disp(max(abs(epipolarConstraintFromE(data.matchEssentialTruth(:,:,1),x1,x2))))
