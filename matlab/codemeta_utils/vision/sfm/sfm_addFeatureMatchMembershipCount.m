%For each feature points, add counts of how many matches it is part of
%function data=sfm_addStructureMatchMembershipCount(data)
%Requires fields feature and match. Adds fields
%data.feature().matchMembershipCount with counts of matches for each point
%and data.feature().matchMembershipMaxCount with the maximum of these.
function data=sfm_addFeatureMatchMembershipCount(data,varargin)
matchMemberName='matchFiltered';
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'member'
            ivarargin=ivarargin+1;
            matchMemberName=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NFeatures=length(data.feature);
NMatches=length(data.(matchMemberName));
for iFeature=1:NFeatures
    NPoints=size(data.feature(iFeature).location,2);
    data.feature(iFeature).matchMembershipCount=zeros(1,NPoints);
end

%for each match, get the two corresponding feature indeces and increase
%counters
for iMatch=1:NMatches
    idxImg=data.(matchMemberName)(iMatch).idxImg;
    m=data.(matchMemberName)(iMatch).idxMatch;
    for jImg=1:2
        data.feature(idxImg(jImg)).matchMembershipCount(m(jImg,:))=...
            data.feature(idxImg(jImg)).matchMembershipCount(m(jImg,:))+1;
    end
end

%store maximum
for iFeature=1:NFeatures
    data.feature(iFeature).matchMembershipMaxCount=max(data.feature(iFeature).matchMembershipCount);
end
