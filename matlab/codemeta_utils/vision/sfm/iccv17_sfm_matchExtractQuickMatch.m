function data=iccv17_sfm_matchExtractQuickMatch(data,optsQuickMatch,varargin)
methodQuickMatch='standard';

%parse optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'nbest'
            methodQuickMatch='nbest';
            ivarargin=ivarargin+1;
            NBest=varargin{ivarargin};
        case 'hierarchical'
            methodQuickMatch='hierarchical';
            ivarargin=ivarargin+1;
            optsQuickMatchHierarchical=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

matchMemberName='matchQuickMatch';
NImages=sfm_getImageNumber(data);

%add default options (notice that they can be overridden by input)
optsQuickMatch = ['ratioInterCluster' 0.75 optsQuickMatch];

%get data from all images to allow joint matching
dataJointMatch=sfm_featureExtractJointMatchData(data,'NImages',NImages);

descriptors=double(dataJointMatch.descriptor);
membershipPrior=dataJointMatch.membershipPrior;

%do the joint matching with QuickMatch
switch methodQuickMatch
    case 'standard'
        D=sqrt(euclideanDistMatrix(descriptors));
        membershipMatch=quickshift_matching(D,membershipPrior,optsQuickMatch{:});
    case 'nbest'
        [D,DIndicator]=matchingEuclideanDistMatrix(descriptors,...
            'membershipPrior',membershipPrior,...
            'NBest',NBest);
        D=sqrt(D);
        membershipMatch=quickshift_matching(D,membershipPrior,...
            'relationIndicator',DIndicator,optsQuickMatch{:});
    case 'hierarchical'
        membershipMatch=quickshift_matchingHierarchical(...
            descriptors,membershipPrior,...
            optsQuickMatchHierarchical{:},'optsMatching',optsQuickMatch);
    otherwise
        error('Variant of QuickMatch not recognized')
end
    
%transform multimach data into pairwise match structure
match=iccv17_clusterMembership2pairwiseMatches(membershipPrior,membershipMatch);

%store result back in the sfm data structure
data.(matchMemberName)=match;
