function EWalk=walksRemoveRedundant(EWalk,varargin)

flagRemoveReverse=true;
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'flagremovereverse'
            ivarargin=ivarargin+1;
            flagRemoveReverse=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

iWalk=1;
while iWalk<size(EWalk,1)
    flagRepetition=findRepetitions(EWalk(iWalk,:),EWalk(iWalk+1:end,:));
    if flagRemoveReverse
        flagRepetition=flagRepetition | ...
            findRepetitions(fliplr(EWalk(iWalk,:)),EWalk(iWalk+1:end,:));
    end
    EWalk([false(iWalk,1); flagRepetition],:)=[];
    iWalk=iWalk+1;
end

function flagRepetition=findRepetitions(thisWalk,EWalk)
[NWalksNext,K]=size(EWalk);
flagRepetition=false(NWalksNext,1);
for ik=1:K
    flagRepetition=flagRepetition ...
        | all(EWalk==repmat(thisWalk,NWalksNext,1),2);
    thisWalk=circshift(thisWalk,[0 1]);
end
