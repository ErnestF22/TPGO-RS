function sfm_displayMatchResiduals(data,varargin)
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
NMatch=length(iMatchList);
NCols=ceil(sqrt(NMatch));
NRows=ceil(NMatch/NCols);
for iiMatch=1:NMatch
    iMatch=iMatchList(iiMatch);
    subplot(NCols,NRows,iiMatch)
    e=data.(matchMemberName)(iMatch).residuals;
    th=data.(matchMemberName)(iMatch).threshold;
    thAuto=sfm_rawThresholdEstimate(e);
    cumDistPerc(e);
    hold on
    plot([th th],[0 100],'r')
    plot([thAuto thAuto],[0 100],'r:')
    hold off
end
