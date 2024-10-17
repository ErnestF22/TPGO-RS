%Transform multimatch clusters to pairwise matches
%function match=iccv17_clusterMembership2pairwiseMatches(membershipPrior,membeshipMatch,NImages)
%Output arguments
%   pMatch  a [NImages*(NImages-1)/2 x 1] array of structures with fields
%       idxImg  [2 x 1] array with indeces of images for that pairwise match
%       idxMatch [2 x m] array with indeces of features in each image that
%           correspond to each other
function match=iccv17_clusterMembership2pairwiseMatches(membershipPrior,membershipMatch)
%number of images
NImages=length(unique(membershipPrior));

%prepare data structure for pairwise matches
NMatchPairs=NImages*(NImages-1)/2;
match=repmat(struct('idxImg',[],'idxMatch',[]),1,NMatchPairs);
cntMatch=1;

%fill in data structure
for iImg = 1:NImages-1
    for jImg = iImg+1:NImages
        %indeces of images
        match(cntMatch).idxImg=[iImg;jImg];
        
        %get membershipMatch for each himage
        membershipMatchIImg=membershipMatch(membershipPrior==iImg);
        membershipMatchJImg=membershipMatch(membershipPrior==jImg);

        %extract match information and put it in the structure
        [~,idxIImg,idxJImg]=intersect(membershipMatchIImg,membershipMatchJImg);
        match(cntMatch).idxMatch = [idxIImg idxJImg]';
        
        %update counter
        cntMatch=cntMatch+1;
    end
end

