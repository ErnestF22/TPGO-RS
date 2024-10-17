function sfm_displayMatch(data,varargin)
optsDisplayMatch={};
flagGivenMatchList=false;
matchMemberName='match';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'matchlist'
            ivarargin=ivarargin+1;
            flagGivenMatchList=true;
            iMatchList=varargin{ivarargin};
        case 'member'
            ivarargin=ivarargin+1;
            matchMemberName=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if ~flagGivenMatchList
    iMatchList=1:length(data.(matchMemberName));
end

matchIdxImg=sfm_getMatchIdxImg(data,matchMemberName);

NMatch=length(iMatchList);
NCols=ceil(sqrt(NMatch));
NRows=ceil(NMatch/NCols);
if NMatch>1
    optsDisplayMatch=[optsDisplayMatch 'noSlider'];
end
for iiMatch=1:NMatch
    iMatch=iMatchList(iiMatch);
    idxImg=matchIdxImg(:,iMatch);
    img1=sfm_getImageById(data,idxImg(1));
    img2=sfm_getImageById(data,idxImg(2));
    [x1,x2]=sfm_getFeatureLocationFromMatchId(data,iMatch,'member',matchMemberName);
    subplot(NCols,NRows,iiMatch)
    sfm_rawDisplayMatch(img1,img2,x1,x2,optsDisplayMatch{:})
end
