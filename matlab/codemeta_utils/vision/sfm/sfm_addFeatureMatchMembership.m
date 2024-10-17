%Add, for each feature point, list of what other feature it correspond
%function data=sfm_addFeatureMatchMembership(data)
%Requires fields feature and match. Adds these sub-fields to feature()
%   matchMembershipCount, matchMembershipMaxCount 
%       see sfm_addFeatureMatchMembershipCount
%   matchMembershipListFeature  
%       indeces of the images where this feature also appears
%   matchMembershipListIdxFeature
%       indeces under which this feature appears in the other images
function data=sfm_addFeatureMatchMembership(data)
memberName='matchFiltered';

data=sfm_addFeatureMatchMembershipCount(data,'member',memberName);

NFeatures=length(data.feature);
NMatches=length(data.(memberName));

%allocate space for membership lists
for iFeature=1:NFeatures
    NPoints=size(data.feature(iFeature).location,2);
    NList=data.feature(iFeature).matchMembershipMaxCount;
    data.feature(iFeature).matchMembershipListFeature=zeros(NList,NPoints);
    data.feature(iFeature).matchMembershipListIdxFeature=zeros(NList,NPoints);
end

%loop over all image pairs
for iMatch=1:NMatches
    idxImg=data.(memberName)(iMatch).idxImg;
    m=data.(memberName)(iMatch).idxMatch;
    Nm=size(m,2);
    %loop source image over first and second image
    for jImg=1:2
        jpImg=mod(jImg,2)+1;   %maps 1=>2 and 2=>1
        %get indexes of source and dest images
        idxImgSource=idxImg(jImg);
        idxImgDest=idxImg(jpImg);
        %loop over feature matches
        for im=1:Nm
            idxFeatureSource=m(jImg,im);
            idxFeatureDest=m(jpImg,im);
            %go to feature in source image, and add destination image, and
            %index in the destination image to the membership list
            data=addFeatureToMembershipList(data,idxImgSource,idxImgDest,idxFeatureSource,idxFeatureDest);
        end
    end
end

function data=addFeatureToMembershipList(data,idxImgSource,idxImgDest,idxFeatureSource,idxFeatureDest)
%find next empty place in the list
idxEmpty=find(data.feature(idxImgSource).matchMembershipListFeature(:,idxFeatureSource)==0,1,'first');

data.feature(idxImgSource).matchMembershipListFeature(idxEmpty,idxFeatureSource)=idxImgDest;
data.feature(idxImgSource).matchMembershipListIdxFeature(idxEmpty,idxFeatureSource)=idxFeatureDest;
